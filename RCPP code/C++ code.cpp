#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

// Define a structure for the phylogenetic tree
struct PhyloTree {
    std::vector<std::string> tip_labels;
    std::vector<std::string> node_labels;
    std::vector<std::vector<int>> edges;
    std::vector<double> edge_lengths;
    int Nnode;
};

// Function to calculate total abundance descending from each branch of a phylo object
std::map<int, double> abundance_phylo(const PhyloTree& tree, const std::vector<std::pair<std::string, double>>& abundance_data) {
    std::vector<std::string> node_labels = tree.tip_labels;
    node_labels.insert(node_labels.end(), tree.node_labels.begin(), tree.node_labels.end());

    std::map<int, double> total_abundance;

    std::function<double(int)> calculate_abundance = [&](int node) {
        double node_abundance = 0;
        for (const auto& data : abundance_data) {
            if (data.first == node_labels[node - 1]) {  // Adjusting for 1-based indexing
                node_abundance = data.second;
                break;
            }
        }

        for (const auto& edge : tree.edges) {
            if (edge[0] == node) {
                double child_abundance = 0;
                if (std::find(tree.tip_labels.begin(), tree.tip_labels.end(), node_labels[edge[1] - 1]) == tree.tip_labels.end()) {
                    child_abundance = calculate_abundance(edge[1]);
                } else {
                    for (const auto& data : abundance_data) {
                        if (data.first == node_labels[edge[1] - 1]) {
                            child_abundance = data.second;
                            total_abundance[edge[1]] = data.second;
                            break;
                        }
                    }
                }
                node_abundance += child_abundance;
            }
        }

        total_abundance[node] = node_abundance;
        return node_abundance;
    };

    calculate_abundance(tree.edges[0][0]);
    return total_abundance;
}

// Function to calculate all ancestor nodes of a given node, including itself
std::vector<int> getAllAncestors(const PhyloTree& tree, int node) {
    int root_node = tree.tip_labels.size() + 1;
    std::vector<int> ancestors = {node};
    int anc_node = root_node + 1;

    if (node != root_node) {
        while (anc_node > root_node) {
            auto it = std::find_if(tree.edges.begin(), tree.edges.end(), [&](const std::vector<int>& edge) {
                return edge[1] == node;
            });
            if (it != tree.edges.end()) {
                anc_node = (*it)[0];
                ancestors.push_back(anc_node);
                node = anc_node;
            } else {
                break;
            }
        }
    }

    return ancestors;
}

// Function to calculate all descendant branches of a node, recording distance from node or "x"
std::vector<std::vector<double>> get_descendant_branches(const PhyloTree& tree, int node) {
    bool root = false;
    std::vector<int> branches;
    if (node == tree.tip_labels.size() + 1) {
        root = true;
        for (size_t i = 0; i < tree.edges.size(); ++i) {
            if (tree.edges[i][0] == node) {
                branches.push_back(i);
            }
        }
    } else {
        for (size_t i = 0; i < tree.edges.size(); ++i) {
            if (tree.edges[i][1] == node) {
                branches.push_back(i);
                break;
            }
        }
    }

    std::vector<std::vector<double>> branch_info;

    for (const auto& branch : branches) {
        std::vector<double> branch_data = {
            static_cast<double>(tree.edges[branch][0]),
            static_cast<double>(tree.edges[branch][1]),
            tree.edge_lengths[branch],
            0.0  // Placeholder for distance
        };
        branch_info.push_back(branch_data);
    }

    for (auto& info : branch_info) {
        if (!root) {
            info[3] = std::find_if(tree.edges.begin(), tree.edges.end(), [&](const std::vector<int>& edge) {
                return edge[1] == node;
            }) != tree.edges.end() ? 1.0 : 0.0;
        }
    }

    return branch_info;
}

// Function to calculate T_i and S_i for a given node i in a phylogenetic tree
std::pair<double, std::vector<std::pair<double, double>>> compute_T_i_S_i(const PhyloTree& tree, int node, const std::map<int, double>& abundances) {
    std::vector<std::vector<double>> df_descen_info;

    for (const auto& edge : tree.edges) {
        if (edge[0] == node) {
            df_descen_info.push_back({static_cast<double>(edge[0]), static_cast<double>(edge[1]), tree.edge_lengths[&edge - &tree.edges[0]]});
        }
    }

    double T_i = 0;
    double prev_x = 0;
    std::vector<std::pair<double, double>> df_list;

    while (!df_descen_info.empty()) {
        double branch_length = df_descen_info.front()[2];
        double abundance_sum = 0;

        for (const auto& info : df_descen_info) {
            abundance_sum += abundances.at(static_cast<int>(info[1]));
        }

        df_list.push_back({abundance_sum, branch_length});
        T_i += abundance_sum * (branch_length - prev_x);
        prev_x = branch_length;

        df_descen_info.erase(std::remove_if(df_descen_info.begin(), df_descen_info.end(), [&](const std::vector<double>& info) {
            return info[2] == branch_length;
        }), df_descen_info.end());
    }

    return {T_i, df_list};
}

// Function to calculate S_i_a for a given node i and ancestor a in a phylogenetic tree
std::pair<std::vector<std::pair<double, double>>, std::vector<double>> calculate_S_i_a(const PhyloTree& tree, int node, const std::map<int, double>& abundances, int curr_ancestor, double h, double l_i) {
    std::vector<std::vector<double>> df_node_info;

    if (node != curr_ancestor) {
        df_node_info = get_descendant_branches(tree, curr_ancestor);
        df_node_info.erase(std::remove_if(df_node_info.begin(), df_node_info.end(), [&](const std::vector<double>& info) {
            return info[3] <= h;
        }), df_node_info.end());

        for (auto& info : df_node_info) {
            if (info[3] - info[2] < h) {
                info[2] = info[3] - h;
            }
            info[3] -= h;
        }
    } else {
        df_node_info = get_descendant_branches(tree, node);
    }

    std::vector<std::pair<double, double>> df_list;
    std::vector<double> abund_list;
    double x = 0;

    while (x != l_i) {
        std::vector<int> index_curr_branches;

        for (size_t i = 0; i < df_node_info.size(); ++i) {
            if (df_node_info[i][3] == df_node_info[i][2]) {
                index_curr_branches.push_back(i);
            }
        }

        if (!df_node_info.empty()) {
            double abundance_sum = 0;
            for (const auto& idx : index_curr_branches) {
                abundance_sum += abundances.at(static_cast<int>(df_node_info[idx][1]));
            }

            df_list.push_back({abundance_sum, *std::min_element(df_node_info.begin(), df_node_info.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
                return a[3] < b[3];
            })[3] + x});
            abund_list.push_back(abundance_sum);

            double prev_x = *std::min_element(df_node_info.begin(), df_node_info.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
                return a[3] < b[3];
            })[3];
            x += prev_x;

            for (auto& info : df_node_info) {
                info[3] -= prev_x;
                if (info[3] == 0) {
                    df_node_info.erase(std::remove(df_node_info.begin(), df_node_info.end(), info), df_node_info.end());
                }
            }
        } else {
            break;
        }
    }

    return {df_list, abund_list};
}

// Function to calculate E/J/M or all given i and a
std::vector<double> calculate_EJM_i_a(const PhyloTree& tree, int node, const std::map<int, double>& abundances, int curr_ancestor, double h, double l_i, char index_letter, int q, bool individual) {
    auto S_i_a_res = calculate_S_i_a(tree, node, abundances, curr_ancestor, h, l_i);
    auto df_S_i_a = S_i_a_res.first;
    auto abund_list = S_i_a_res.second;

    std::vector<double> df_E_i_a, df_J_i_a, df_M_i_a;

    for (size_t k = 0; k < abund_list.size(); ++k) {
        double S = df_S_i_a[k].first;
        double m = std::log(abund_list.size());
        double e = 0, j = 0;

        for (const auto& abund : abund_list) {
            e += -(abund / S) * std::log(abund / S);
        }

        if (abund_list.size() > 1) {
            for (const auto& abund : abund_list) {
                j += -(abund / S) * std::log(abund / S) / std::log(abund_list.size());
            }
        } else {
            j = 1;
        }

        df_M_i_a.push_back(m);
        df_E_i_a.push_back(e);
        df_J_i_a.push_back(j);
    }

    if (individual) {
        return {df_E_i_a[0], df_J_i_a[0], df_M_i_a[0]};
    } else {
        return {df_E_i_a[0], df_J_i_a[0], df_M_i_a[0]};
    }
}

// Function to calculate the integral for each ancestor
double calculate_integral(const PhyloTree& tree, int node, int curr_ancestor, const std::vector<std::pair<double, double>>& S_i, const std::map<int, double>& abundances, char index_letter, int q, bool individual) {
    double l_i = S_i.back().second;
    if (l_i == 0) return 0;

    double h = 0;
    int temp_node = node;
    while (temp_node != curr_ancestor) {
        for (const auto& edge : tree.edges) {
            if (edge[1] == temp_node) {
                h += tree.edge_lengths[&edge - &tree.edges[0]];
                temp_node = edge[0];
                break;
            }
        }
    }

    double int_2 = 0;
    double prev_x = h;
    size_t current_row = 0;
    while (current_row < S_i.size() && S_i[current_row].second > h) {
        double x_s = S_i[current_row].second;
        double value_s = S_i[current_row].first;
        if (x_s > l_i) {
            int_2 += value_s * (l_i - prev_x);
            break;
        } else {
            int_2 += value_s * (x_s - prev_x);
            prev_x = x_s;
            ++current_row;
        }
    }

    std::function<double(const std::vector<double>&)> integral = [&](const std::vector<double>& df) {
        double prev_x = 0;
        size_t current_row_s = 0, current_row_i = 0;
        double int_1 = 0;

        while (current_row_s < S_i.size() && current_row_i < df.size()) {
            double x_s = S_i[current_row_s].second;
            double x_i = df[current_row_i];
            double value_s = S_i[current_row_s].first;
            double value_i = df[current_row_i];

            if (x_s < x_i) {
                int_1 += value_s * value_i * (x_s - prev_x);
                prev_x = x_s;
                ++current_row_s;
            } else if (x_s == x_i) {
                int_1 += value_s * value_i * (x_s - prev_x);
                prev_x = x_s;
                ++current_row_s;
                ++current_row_i;
            } else {
                int_1 += value_s * value_i * (x_i - prev_x);
                prev_x = x_i;
                ++current_row_i;
            }
        }

        return int_1;
    };

    if (int_2 != 0) {
        auto Sum_i_a = calculate_EJM_i_a(tree, node, abundances, curr_ancestor, h, l_i, index_letter, q, individual);
        return integral(Sum_i_a) * int_2;
    } else {
        return 0;
    }
}

// Function to calculate E/J/M for a given node
std::vector<double> calculate_EJM_i(const PhyloTree& tree, int node, const std::map<int, double>& abundances, char index_letter, int q, bool individual) {
    auto S_i_all = compute_T_i_S_i(tree, node, abundances);
    double T_i = S_i_all.first;
    auto S_i = S_i_all.second;

    auto ancestors = getAllAncestors(tree, node);
    double ejm_i = 0;
    for (const auto& ancestor : ancestors) {
        ejm_i += calculate_integral(tree, node, ancestor, S_i, abundances, index_letter, q, individual);
    }

    if (T_i != 0) {
        return {(1 / T_i) * ejm_i};
    } else {
        return {0};
    }
}

// Function to calculate D/J indices for a given node
std::vector<double> calculate_DJ_i(const PhyloTree& tree, int node, const std::map<int, double>& abundances, char index_letter, int q, bool individual) {
    auto ancestors = getAllAncestors(tree, node);
    auto S_i_all = compute_T_i_S_i(tree, node, abundances);
    double T_i = S_i_all.first;
    auto S_i = S_i_all.second;

    if (T_i != 0) {
        double dj_i = 0;
        for (const auto& ancestor : ancestors) {
            dj_i += calculate_integral(tree, node, ancestor, S_i, abundances, index_letter, q, individual);
        }

        return {(1 / T_i) * dj_i};
    } else {
        return {0};
    }
}

// Function to calculate longitudinal or star mean for D for q=0 or 1, or J for q=1, or all for one mean type
std::vector<double> long_star(const PhyloTree& tree, const std::vector<std::pair<std::string, double>>& node_abundances, const std::string& mean_type, char index_letter, int q, bool individual) {
    auto abundances = abundance_phylo(tree, node_abundances);
    int node = tree.edges[0][0];

    auto S_i_a_res = calculate_S_i_a(tree, node, abundances, node, 0, std::accumulate(tree.edge_lengths.begin(), tree.edge_lengths.end(), 0.0));
    auto df_S_i_a = S_i_a_res.first;
    auto abund_list = S_i_a_res.second;

    std::vector<double> df_E_i_a, df_J_i_a, df_M_i_a;

    for (size_t k = 0; k < abund_list.size(); ++k) {
        double S = df_S_i_a[k].first;
        double m = std::log(abund_list.size());
        double e = 0, j = 0;

        for (const auto& abund : abund_list) {
            e += -(abund / S) * std::log(abund / S);
        }

        if (abund_list.size() > 1) {
            for (const auto& abund : abund_list) {
                j += -(abund / S) * std::log(abund / S) / std::log(abund_list.size());
            }
        } else {
            j = 1;
        }

        df_M_i_a.push_back(m);
        df_E_i_a.push_back(e);
        df_J_i_a.push_back(j);
    }

    if (individual) {
        return {df_E_i_a[0], df_J_i_a[0], df_M_i_a[0]};
    } else {
        return {df_E_i_a[0], df_J_i_a[0], df_M_i_a[0]};
    }
}

// Function to calculate distance
double distance(PhyloTree& tree, int top_node, int bottom_node) {
    int curr_node = bottom_node; // Select bottom node as starting node
    double dist = 0.0; // Initialize distance sum

    while (curr_node != top_node) { // Do until reach desired top node
        // Each node only has one parent branch,
        // select that branch, add the distance between node and parent,
        // move to parent node and repeat
        auto it = std::find_if(tree.edges.begin(), tree.edges.end(), [curr_node](const std::pair<int, int>& edge) {
            return edge.second == curr_node;
        });
        if (it != tree.edges.end()) {
            int index = std::distance(tree.edges.begin(), it);
            dist += tree.edge_lengths[index]; // Sum distance
            curr_node = tree.edges[index].first; // Move to parent node
        }
    }
    return dist;
}

// Function to calculate node-wise mean for D for q = 0 or 1, or for J for q = 1, or all
std::unordered_map<std::string, double> node(const std::string& file, bool node_abundances = false, std::string index_letter = "D", int q = 1, bool individual = false) {
    PhyloTree tree;
    // Read and convert file to tree (This function needs to be implemented based on the file format)
    read_convert(file, tree);

    std::transform(index_letter.begin(), index_letter.end(), index_letter.begin(), ::toupper);

    // Check if the tree is linear
    if (tree.tip_label.size() == 1) {
        if (individual) {
            return { {"index", 1.0} };
        } else {
            return { {"D1N", 1.0}, {"J1N", 1.0}, {"D0N", 1.0} };
        }
    } else {
        // Check if tree has abundance data
        if (!tree.node_abundances.empty()) {
            abundances = abundance_phylo(tree); // Calculate branch/node abundances
        } else {
            int num_tips = tree.edge[0].first - 1; // Number of tips
            for (int i = 1; i <= num_tips; ++i) {
                tree.tip_label.push_back(std::to_string(i)); // Assign node labels
            }

            int num_nodes = tree.Nnode; // Number of nodes
            for (int i = num_tips + 1; i <= num_tips + num_nodes; ++i) {
                tree.node_abundances[std::to_string(i)] = 0;
            }

            for (int i = 1; i <= num_tips; ++i) {
                tree.node_abundances[std::to_string(i)] = 1.0 / num_tips;
            }

            abundances = abundance_phylo(tree); // Calculate branch/node abundances
        }

        // Calculate h_bar
        double T = 0.0; // Initialize sum
        for (size_t i = 0; i < tree.edge.size(); ++i) {
            T += tree.node_abundances[std::to_string(tree.edge[i].second)] * tree.edge_length[i];
        }

        // Numbers of every internal node, corresponding to numbers in phylo object
        std::vector<int> nodes(tree.tip_label.size() + 1, tree.tip_label.size() + tree.Nnode);

        // Calculate node-averaged indices
        if (!individual) { // All indices
            std::vector<std::vector<double>> DJ;
            for (int node : nodes) {
                DJ.push_back(calculate_DJ_i(node, tree, abundances, individual, q, "D"));
            }

            double D1N = 0.0, J1N = 0.0, D0N = 0.0;
            for (const auto& dj : DJ) {
                D1N += dj[0];
                J1N += dj[1];
                D0N += dj[2];
            }
            D1N = std::exp(D1N / T);
            J1N /= T;
            D0N = std::exp(D0N / T);

            return { {"D1N", D1N}, {"J1N", J1N}, {"D0N", D0N} };
        } else { // One index
            std::vector<double> DJ;
            for (int node : nodes) {
                DJ.push_back(calculate_DJ_i(node, tree, abundances, individual, q, index_letter));
            }

            double index = 0.0;
            if (index_letter == "J") {
                index = std::accumulate(DJ.begin(), DJ.end(), 0.0) / T;
            } else {
                index = std::exp(std::accumulate(DJ.begin(), DJ.end(), 0.0) / T);
            }

            return { {"index", index} };
        }
    }
}
