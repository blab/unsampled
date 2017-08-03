#import packages
from dendropy import Tree, Node, Edge
import math
import numpy as np
import random as rd
import operator
from scipy.stats import zscore

#function that takes the product of multiple numbers, or returns 1 if the list is empty
def prod(factors):
    return reduce(operator.mul, factors, 1)

#names nodes in a preorder tree traversal order (after tips are named) for clarity during function testing
def name_nodes(tree):
    current_node_name = len(tree.leaf_nodes())+1
    for (index, node) in enumerate(tree.preorder_node_iter()):
        if node.taxon:
            node.label = str(node.taxon.label)
        else:
            node.label = str(current_node_name)
            current_node_name += 1

def name_edges(tree):
    #give each edge a label
    for (index, edge) in enumerate(tree.preorder_edge_iter()):
        edge.label = str(index)

def zipped_sorted_intervals(tree):
    """
    Takes in a tree object and finds information on each interval (marked by time points from root).
    Returns a list of lists. Each element in the final list represents an interval. Within each interval is:
        1. [start time of interval, end time of interval]
        2. list of active lineages within each interval (referenced by the node at its tail)
    """
    #each node and its distance from root
    zipped_dists = zip(tree.nodes(),tree.calc_node_root_distances(return_leaf_distances_only =False))
    sorted_zipped_dists = reversed(sorted(zipped_dists, key=lambda branches: branches[1])) #sort by farthest to nearest from root

    intervals = []
    living_lineages = []
    current_start = tree.max_distance_from_root()

    for (node, distance) in sorted_zipped_dists:
        if current_start == distance: #if there are multiple lineages sampled at the same time, add them to the current interval and move to the next node
            living_lineages.append(node)

        else:
            intervals.append([[current_start, distance], living_lineages]) #add the interval to the set of intervals
            current_start = distance #update the starting distance of the new interval

            children =set(node.child_nodes()) #add the current node and remove children of the node if there are any
            living_lineages = list(children.symmetric_difference(living_lineages))
            living_lineages.append(node)
    return intervals

def zipped_partial_intervals(tree, cut_dist_from_root):
    """
    Like zipped_sorted_intervals but takes a cut_distance_from_root, the time point in the tree
    to treat as time 0 (ignores all of the branches and samples that occur more recently from this time)
    """
    zipped_lineages = zipped_sorted_intervals(tree)

    if cut_dist_from_root > tree.max_distance_from_root():
        included_intervals =[([float(cut_dist_from_root), tree.max_distance_from_root()],[])]
        included_intervals.extend(zipped_lineages)

    else:
        included_intervals = [(interval_endpoints, interval_nodes) for (interval_endpoints, interval_nodes) in zipped_lineages if cut_dist_from_root > interval_endpoints[1]]
        included_intervals[0][0][0] = float(cut_dist_from_root) #if the cut fell in the middle of an interval, change the first included interval to reflect this

    return included_intervals

def conditioned_prob_lineage_coal(tree, popsize, cut_dist_from_root = None):
    """
        Takes in a tree object and a constant population size and returns a list containing
        the probability of coalescence to one lineage for each interval.
        Those probabilities are conditioned by the probability of not coalescing anywhere else downstream in the tree.
    """
    if cut_dist_from_root:
        intervals = zipped_partial_intervals(tree, cut_dist_from_root)
    else:
        intervals = zipped_sorted_intervals(tree)

    interval_length = []
    num_lineages = []

    for (interval_endpoints, nodes) in intervals: #for each interval find the length of the interval and the number of lineages present
        interval_length.append(interval_endpoints[0] - interval_endpoints[1])
        num_lineages.append(len(nodes))

    #continuous time probability of coalescence for each interval (not conditioned on time)
    prob_coalescence_in_interval = [1-math.exp(-(float(lineages)/popsize)*length) if num_lineages is not 0 else 0 for (length, lineages) in zip(interval_length, num_lineages)]
    prob_no_coalescence_in_interval = [1-pcoal for pcoal in prob_coalescence_in_interval]

    coalintervals = []
    #find the probability of coalescing with a single lineage conditioned by position in tree (p no coal earlier * pcoal current interval)
    for (index, interval) in enumerate(prob_no_coalescence_in_interval):
        if num_lineages[index] is not 0:
            coalintervals.append(prod(prob_no_coalescence_in_interval[:index])*prob_coalescence_in_interval[index]/num_lineages[index])#the number of lineages
        else:
            coalintervals.append(0)
    return coalintervals

def pcoal_along_edge(tree, popsize, cut_dist_from_root = None):
    """
        Takes in a tree object and constant population size and returns a dictionary containing
        each edge's conditioned probability of coalescence.
        Each probability is a value keyed to the corresponding edge and these probabilities are not cumulative
    """
    edge_prob = {}
    for node in tree.nodes():
        edge_prob[node.edge] = 0.0

    pcoal = conditioned_prob_lineage_coal(tree, popsize, cut_dist_from_root)

    if cut_dist_from_root:
        interval_lineages = zip(*zipped_partial_intervals(tree, cut_dist_from_root))[1]
    else:
        interval_lineages = zip(*zipped_sorted_intervals(tree))[1]

    for (index, interval) in enumerate(interval_lineages): #for the list of nodes in each interval, if it doesn't belong to the dictionary, set the edge connection to it to pcoal (conditioned)
        for node in interval:
            edge_prob[node.edge] += pcoal[index]

    edge_prob[tree.seed_node.edge] = 1-sum(edge_prob.values())
    return edge_prob

def calculate_cumulative_node_prob(tree, popsize):
    """
        Takes in a tree object and constant population size and finds the cumulative expectation of proportion of the
        tips theoretically sampled under each internal node using the coalescent model.
    """
    prob_lineage = pcoal_along_edge(tree, popsize)

    #dictionary to store the cumulative probability of a new sample coalescing (value) under each node (key)
    cumulative_node_prob ={}

    #look at each internal (non-tip) node
    for node in tree.internal_nodes():
        node_prob = 0

        #iterate through the nodes belonging to the subtree rooted at node
        for subtree_node in node.preorder_iter():

            #look at each edge of the node and add its probability to the cumulative node prob
            for edge in subtree_node.child_edge_iter():
                node_prob += prob_lineage[edge]

        cumulative_node_prob[node] = node_prob
    return cumulative_node_prob

def node_zscores(node_prob_dict, tree):
    """
    Takes in the dictionary of cumulative expected sampling proportions (from calculate_cumulative_node_prob) and the tree
    and gives a measure of z-score for each internal node (undersampled are < 0)
    """
    keys, vals = zip(*node_prob_dict.items())
    num_tips = [len(node.leaf_nodes()) for node in keys]
    total_tips = len(tree.leaf_nodes())
    vals = [prob*total_tips-tips for (prob, tips)  in zip(vals, num_tips)]
    return dict(zip(keys, zscore(vals, ddof=1)))

def randomly_prune_tree(tree, proportion_tips_to_drop):
    """
        Takes in a tree and a proportion of tips to drop and returns, in a dict, the full tree, remaining pruned tree
        list of dropped samples, dictionary of distances from each node to the root of the full tree, and the
        root reference, the distance in the full tree where the root of the pruned tree lies
    """
    name_edges(tree)
    name_nodes(tree)
    pruned_tree = tree.clone()
    name_edges(pruned_tree)
    name_nodes(pruned_tree)

    #take the distance from the original root (to give time references of when to look for a missed sampling event)
    full_dist_from_root = {}
    for nd in pruned_tree.nodes():
        full_dist_from_root[nd] = nd.distance_from_root()

    #drop some proportion of samples
    dropped_samples = rd.sample(pruned_tree.leaf_nodes(), int(round(proportion_tips_to_drop*len(tree.leaf_nodes())))) #actual samples to be dropped
    for nd in dropped_samples:
        pruned_tree.prune_subtree(nd) #prune the tree by dropping tips in dropped_samples

    root_ref = full_dist_from_root.get(pruned_tree.seed_node)

    return {'full_tree': tree, 'pruned_tree': pruned_tree, 'dropped_samples': dropped_samples, 'node_root_dist_full_tree': full_dist_from_root, 'root_reference': root_ref}

def theoretical_attachment_prob_matrix(tree, pruned_tree, root_ref, popsize, full_dist_from_root, dropped_samples):
    """
        Uses the coalescent framework to find the theoretical probability of each dropped sample attaching to each edge in the
        remaining pruned tree. It requires information about both the full and pruned tree and returns a numpy array
        where each row corresponds to a dropped sample and each column represents an edge in the pruned tree. The cells
        are the probabilities that the node attached at that edge and each row sums to 1.
    """
    unsampled_theoretical_probs =[]

    for node in dropped_samples:
        current_node_probs = []
        samp_dist = full_dist_from_root.get(node) - root_ref #makes sure if the root of the original tree was removed, this distance reflects the change
        if samp_dist >= 0: #if the sample is more recent than the pruned tree's root
            probs = pcoal_along_edge(pruned_tree, popsize, samp_dist) #the edge prob coals at the sample's distance from the root

            for edge in pruned_tree.preorder_edge_iter(): #in edge order, add each probability to the current node's probability list
                current_node_probs.append(probs.get(edge))

            assert (abs(sum(current_node_probs)-1.0) < .000001), "The sum of the probabilities for this dropped sample is not 1"
            unsampled_theoretical_probs.append(current_node_probs) #give the theoretical probability matrix the probabilities of the current node for each edge

        else: #if the sample predates the pruned tree's root, it must attach at the edge above the root with probability 1
            current_node_probs = [1.0]
            current_node_probs.extend([0.0]*(len(pruned_tree.edges()) - 1))
            unsampled_theoretical_probs.append(current_node_probs)

    return (np.array(unsampled_theoretical_probs)) #matrix: columns are edges, rows are known unsampled tips, values are probabilities

def observed_attachment_prob_matrix(tree, pruned_tree, dropped_samples):
    """
        Takes in the tree, pruned tree, and dropped samples and finds the edge each dropped sample attaches to in the pruned tree
    """
    pruned_tree_node_labels = [node.label for node in pruned_tree.preorder_node_iter()] #all node labels in the pruned tree
    pruned_tree_edge_labels = [edge.label for edge in pruned_tree.preorder_edge_iter()] #all edge labels in the pruned tree
    node_matrix = []

    for sample in dropped_samples: #look at each dropped sample
        current_sample = [0.0]*len(pruned_tree_edge_labels)
        node = next(n for n in tree.nodes() if n.label == sample.label) #find the corresponding node in the original tree
        current_parent = node.parent_node #parent in the full tree
        current_sibling = node.sibling_nodes()[0]
        found_attachment_edge = False

        while not found_attachment_edge: # looking for where it attaches in the pruned tree
            if current_parent.label in pruned_tree_node_labels: #if its parent is in the pruned tree, that is where the sample connects
                edge_to_attach_to = current_parent.edge.label
                found_attachment_edge = True

            else:
                for nd in current_sibling.preorder_iter(): #otherwise look at the sibling node and all of its decendents to see if they are in the pruned tree
                    if nd.label in pruned_tree_node_labels:
                        edge_to_attach_to = nd.edge.label
                        found_attachment_edge = True
                        break
                else: #the sibling, its children, and the parent are not in the pruned tree, move up to the parent's parent and try again
                    current_sibling = current_parent.sibling_nodes()[0]
                    current_parent = current_parent.parent_node

        current_sample[pruned_tree_edge_labels.index(edge_to_attach_to)] = 1.0 #input a 1 at the edge the dropped sample attaches
        assert (sum(current_sample) == 1), "the sum of the probabilities for this dropped sample is not 1"
        node_matrix.append(current_sample)
    return np.array(node_matrix)

def binary_expectations(expectations):
    """
        Takes in the expectation matrix and for each dropped sample assigns a 1 to the most likely edge of attachment and
        0 to all other edges in the pruned tree
    """
    bin_exp = []
    for row in expectations:
        bin_exp.append([1 if x == max(row) else 0 for x in row])
    return np.array(bin_exp)

def y_dist_dict(tree):
    order=[x for x in tree.leaf_node_iter()]## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
    name_order=[x.label for x in order] #tip labels
    skips=[1 for x in order]
    node_positions = {}

    storePlotted=0
    drawn=[] ## drawn keeps track of what's been drawn
    while len(drawn)!=len(tree.nodes()): # keep drawing the tree until everything is drawn
         for node in [x for x in tree.nodes() if x.label not in drawn]: ## iterate through objects that have not been drawn
            if node.is_leaf(): ## if leaf - get position of leaf, draw branch connecting tip to parent node
                y_idx=name_order.index(node.label) ## y position of leaf is given by the order in which tips were visited during the traversal
                y=sum(skips[:y_idx]) ## sum across skips to find y position

                node_positions[node] = y # y coordinates
                drawn.append(node.label) ## remember that this objects has been drawn

            if node.is_internal(): ## if parent is non-root node and y positions of all its children are known
                if len([node_positions.get(c_nd) for c_nd in node.child_nodes() if node_positions.has_key(c_nd)])==len(node.child_nodes()):
                    children_y_coords=[node_positions.get(c_nd) for c_nd in node.child_nodes() if node_positions.has_key(c_nd)] ## get all existing y coordinates of the node
                    y=sum(children_y_coords)/float(len(children_y_coords)) ## internal branch is in the middle of the vertical bar

                    node_positions[node]=y
                    drawn.append(node.label) ## remember that this objects has been drawn

    assert len(drawn)>storePlotted,'Got stuck trying to find y positions of objects'
    storePlotted=len(drawn)
    return node_positions

if __name__ == '__main__':
    main()
