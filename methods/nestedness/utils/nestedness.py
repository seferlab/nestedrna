import numpy as np
from scipy.special import comb
import scipy.special
import pandas as pd


def degree_of_node(com_i, raw_df_leveled, is_degree_by_domains):
  if is_degree_by_domains == False:
    subset = list(set(raw_df_leveled[raw_df_leveled.values == com_i].values.flatten()))
    sub_subset = list(filter(lambda x: x != com_i, subset))
    neighbors_list = set(filter(lambda x: x.startswith("L"), sub_subset))
    k_i = len(neighbors_list)
  elif is_degree_by_domains == True:
    neighbors_list = set(raw_df_leveled[raw_df_leveled.values == com_i].index.values)
    k_i = len(neighbors_list)
  return k_i, neighbors_list


def find_neighbors(k_i, network_fullname_df):
  source_i_sublist = []
  target_i_sublist = []
  for source in network_fullname_df.source.values:
    if source.startswith(k_i + "_"):
      targets = network_fullname_df[network_fullname_df.source == source].target.values
      target_i_sublist = np.append(target_i_sublist, np.squeeze(targets))
  for target in network_fullname_df.target.values:
    if target.startswith(k_i + "_"):
      sources = network_fullname_df[network_fullname_df.target == target].source.values
      source_i_sublist = np.append(source_i_sublist, np.squeeze(sources))
  source_i_filtered_set = set()
  target_i_filtered_set = set()

  # remove branch info
  for word in target_i_sublist:
    indx = word.find("_")
    community = word[0:indx]
    target_i_filtered_set.add(community)
  for word in source_i_sublist:
    indx = word.find("_")
    community = word[0:indx]
    source_i_filtered_set.add(community)
  neig_i = source_i_filtered_set.union(target_i_filtered_set)
  return neig_i


def find_common_neighbors(n_i, n_j):
  neig_i = find_neighbors(n_i)
  k_i = len(neig_i)
  neig_j = find_neighbors(n_j)
  k_j = len(neig_j)
  nodes_shared = neig_i.intersection(neig_j)
  return nodes_shared, k_i, k_j, neig_i, neig_j


def binom(n, k):
  return scipy.special.comb(n, k, exact=True)


def prob_ij(k_i, k_j, n_tot_num_nodes):
  k_degree_min = np.min([k_i, k_j])
  vals = []
  for k in np.arange(1, k_degree_min + 1, 1):
    val_num = binom(n_tot_num_nodes, k) * binom(n_tot_num_nodes - k, k_j - k) * binom(n_tot_num_nodes - k_j, k_i - k)
    val_den = binom(n_tot_num_nodes, k_j) * binom(n_tot_num_nodes, k_i)
    val = k * (val_num / val_den)
    vals = np.append(vals, val)
  return k_degree_min, np.sum(vals)


def p_val(nodes_shared, n_tot_num_nodes, k_i, k_j):
    k_degree_min = np.min([k_i, k_j])
    p_less = 0
    p_exact = 0
    p_greater = 0
    for k in np.arange(0, k_degree_min + 1, 1):
        val_num = binom(n_tot_num_nodes, k) * binom(n_tot_num_nodes - k, k_j - k) * binom(n_tot_num_nodes - k_j, k_i - k)
        val_den = binom(n_tot_num_nodes, k_j) * binom(n_tot_num_nodes, k_i)
        if k < nodes_shared:
          p_less += (val_num / val_den)
        if k == nodes_shared:
          p_exact += (val_num / val_den)
        if k > nodes_shared:
          p_greater += (val_num / val_den)
    return p_less+p_exact, p_exact, p_greater+p_exact


def w_ij(nodes_shared, prob, k_degree_min):
  return (nodes_shared - prob) / k_degree_min


def omega_ij(k_degree_min, prob, shared_nodes, k_i, k_j, n_tot_num_nodes):
  if shared_nodes > prob:
    return (k_degree_min - prob) / k_degree_min
  elif shared_nodes < prob:
    if (k_i + k_j - n_tot_num_nodes) < 0:
      return (prob) / k_degree_min
    else:
      #print(f"expected = {prob}; ki+kj-n = {(k_i + k_j - n_tot_num_nodes)}")
      return (prob - (k_i + k_j - n_tot_num_nodes)) / k_degree_min
  elif shared_nodes == prob:
    return 1


def n_ij(w_coef, omega_coef):
  return w_coef / omega_coef