import argparse
from ete3 import Tree

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

def main():
    parser = argparse.ArgumentParser(description="Split a large phylogenetic tree into subtrees based on a maximum number of tips.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input tree")
    parser.add_argument('-o', '--output', type=str, required=True, help="Directory to save the subtrees' tip names")
    parser.add_argument('-m', '--max-tips', type=int, required=True, help="Maximum number of tips")

    args = parser.parse_args()

    tree = Tree(args.input)
    
    subtrees = split_tree(tree, args.max_tips)
    
    total_subtrees = len(subtrees)
    for i, subtree in enumerate(subtrees):
        tip_names = [leaf.name for leaf in subtree.get_leaves()]
        
        output_file = f"{args.output}/subset_{i+1}_of_{total_subtrees}.txt"
        
        with open(output_file, 'w') as f:
            for tip in tip_names:
                f.write(f"{tip}\n")
        
        print(f"Subset {i+1} of {total_subtrees} saved to {output_file}")

if __name__ == "__main__":
    main()
