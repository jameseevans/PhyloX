import argparse
from ete3 import Tree
from Bio import SeqIO
import os

def split_tree(tree, max_tips):
    """
    Recursively split a tree into monophyletic subtrees.
    """
    def split_node(node):
        if len(node) > max_tips:
            subtrees = []
            for child in node.children:
                if len(child) > 1:
                    subtrees.extend(split_node(child))
            return subtrees
        else:
            return [node]
    return split_node(tree)

def preserve_root_branch_length(subtree, original_branch_length):
    """
    Reattach the original branch length to the root of the subtree.
    """
    subtree.dist = original_branch_length
    return subtree

# Custom function to write trees while keeping full branch lengths to avoid rounding with ETE3
def write_newick_with_precision(tree, filename, precision=15):
    """
    Write tree to Newick format with specified precision for branch lengths.
    """
    def format_node(node):
        if node.children:
            children_str = ",".join(format_node(child) for child in node.children)
            node_str = f"({children_str})"
        else:
            node_str = ""
        
        if hasattr(node, 'name') and node.name:
            node_str += node.name
        
        if hasattr(node, 'dist') and node.dist is not None:
            node_str += f":{node.dist:.{precision}g}"
        
        return node_str
    
    newick_str = format_node(tree) + ";"
    
    with open(filename, 'w') as f:
        f.write(newick_str)

def filter_fasta(input_fasta, output_fasta, tip_names):
    """
    Filter FASTA file to only sequences in tip_names.
    """
    tip_names = set(tip_names)
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in tip_names:
                SeqIO.write(record, out_f, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Split a phylogenetic tree into monophyletic subtrees.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input tree (Newick format)")
    parser.add_argument('-o', '--output', type=str, required=True, help="Directory to save the subtrees")
    parser.add_argument('-m', '--max-tips', type=int, required=True, help="Maximum number of tips in each subtree")
    parser.add_argument('-a', '--alignment', type=str, help="Optional Fasta alignment to filter for each subtree")
    parser.add_argument('-p', '--precision', type=int, default=15, help="Precision for branch lengths (Sig. dig; default: 15)")
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    tree = Tree(args.input, format=0)
    
    subtrees = split_tree(tree, args.max_tips)
    
    for i, node in enumerate(subtrees):
        tip_names = [leaf.name for leaf in node.get_leaves()]
        
        tip_file = os.path.join(args.output, f"subset_{i+1}.txt")
        with open(tip_file, 'w') as f:
            f.writelines(f"{tip}\n" for tip in tip_names)
        
        subtree_copy = node.copy(method='deepcopy')
        subtree_copy = preserve_root_branch_length(subtree_copy, node.dist)
        
        nwk_file = os.path.join(args.output, f"subset_{i+1}.tre")
        write_newick_with_precision(subtree_copy, nwk_file, args.precision)
        
        if args.alignment:
            fasta_file = os.path.join(args.output, f"subset_{i+1}.fasta")
            filter_fasta(args.alignment, fasta_file, tip_names)
        
        print(f"Subset {i+1} saved: TXT ({tip_file}), NWK ({nwk_file})"
              + (f", FASTA ({fasta_file})" if args.alignment else ""))

if __name__ == "__main__":
    main()
