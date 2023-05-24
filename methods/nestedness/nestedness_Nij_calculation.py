import utils.nestedness as nestedness


def test_metric(nodes_shared, k_i, k_j, n_tot_num_nodes):
    print(f"""
    for a network with total number of domains {n_tot_num_nodes}
    we test the domain overlap between:
    community i with degree {k_i}
    and community j with degree {k_j}.
    We observe that they share {nodes_shared} domains.
    """)
    k_degree_min, prob = nestedness.prob_ij(k_i, k_j, n_tot_num_nodes)
    print("1. Expected number of shared nodes:        ", prob)
    w_coef = nestedness.w_ij(nodes_shared, prob, k_degree_min)
    print("2. calculate w :     ", w_coef)
    omega_coef = nestedness.omega_ij(k_degree_min, prob, nodes_shared, k_i, k_j, n_tot_num_nodes)
    print("3. then Omega:      ", omega_coef)
    n_coef = nestedness.n_ij(w_coef, omega_coef)
    print("4. eventually, Nestedness, N:      ", n_coef)
    p_less, p_exact, p_greater = nestedness.p_val(nodes_shared, n_tot_num_nodes, k_i, k_j)
    print(f"5. Probabilities, less: {p_less}; exact: {p_exact}; greater: {p_greater}")


nodes_shared, k_i, k_j, n_tot_num_nodes = [1, 7, 6, 93]
test_metric(nodes_shared, k_i, k_j, n_tot_num_nodes)
