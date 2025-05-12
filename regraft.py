### Under development, not working for some topologies ###

import argparse
from Bio import Phylo
import os

def write_newick(tree, output_path):
    """Custom function to write Newick with full branch lengths."""
    def get_newick(clade):
        """Recursively generate Newick format for the clade."""
        if clade.is_terminal():
            return f"{clade.name}:{clade.branch_length:.16f}" if clade.branch_length else clade.name
        else:
            subclades = ",".join(get_newick(sub) for sub in clade.clades)
            length = f":{clade.branch_length:.16f}" if clade.branch_length else ""
            return f"({subclades}){length}"
    
    newick_str = get_newick(tree.root) + ";"
    with open(output_path, "w") as file:
        file.write(newick_str)

def graft_subtrees(pruned_tree_path, subtree_dir, output_path):
    pruned_tree = Phylo.read(pruned_tree_path, "newick")
    pruned_tips = {clade.name for clade in pruned_tree.get_terminals()}

    for subtree_file in os.listdir(subtree_dir):
        if subtree_file.endswith(".tre"):
            subtree_path = os.path.join(subtree_dir, subtree_file)
            subtree = Phylo.read(subtree_path, "newick")

            matching_tip = None
            for tip in subtree.get_terminals():
                if tip.name in pruned_tips:
                    matching_tip = tip.name
                    break
                  
            if matching_tip is None:
                continue
            replaced = False
            for parent in pruned_tree.find_clades(order="level"):
                for i, clade in enumerate(parent.clades):
                    if clade.name == matching_tip:
                        parent.clades[i] = subtree.root
                        replaced = True
                        break
                if replaced:
                    break

    write_newick(pruned_tree, output_path)
    print(f"Grafted tree saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Graft subtrees onto a pruned phylogenetic tree.")
    parser.add_argument("-t", "--tree", required=True, help="Path to the pruned tree file (Newick format).")
    parser.add_argument("-s", "--subtrees", required=True, help="Directory containing subtree files (.tre).")
    parser.add_argument("-o", "--output", required=True, help="Output file path for the grafted tree.")
    args = parser.parse_args()

    graft_subtrees(args.tree, args.subtrees, args.output)

if __name__ == "__main__":
    main()
