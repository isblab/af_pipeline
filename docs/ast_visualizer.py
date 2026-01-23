import ast
import graphviz
import json
from ast2json import ast2json

def visualize_ast(code):
    tree = ast.parse(code)
    dot = graphviz.Digraph(comment='Python AST')

    # Recursive function to add nodes and edges
    def add_nodes_edges(node, parent_name=None):
        node_name = str(id(node)) # Unique ID for each node
        node_label = type(node).__name__
        print(node.__dict__)
        dot.node(node_name, node_label)

        if parent_name:
            dot.edge(parent_name, node_name)

        for field, value in ast.iter_fields(node):
            if isinstance(value, list):
                for item in value:
                    if isinstance(item, ast.AST):
                        add_nodes_edges(item, node_name)
            elif isinstance(value, ast.AST):
                add_nodes_edges(value, node_name)

    add_nodes_edges(tree)
    dot.render('ast_visualization', view=True) # Renders and opens the image


def explore_ast(code):

    tree = ast.parse(code)
    for node in ast.walk(tree):
        print(f"Node Type: {type(node).__name__}")
        for field, value in ast.iter_fields(node):
            print(f"  {field}: {value}")
        print("-" * 40)


def ast_to_json(code):
    tree = ast.parse(code)
    ast_json = ast2json(tree)
    print(json.dumps(ast_json, indent=4))

script_path = "/home/omkar/Projects/af_pipeline/af_pipeline/_initialize.py"

with open(script_path, "r") as file:
    code_snippet = file.read()

# Example usage
# code_snippet = """
# def greet(name):
#     print(f"Hello, {name}!")
# """
# visualize_ast(code_snippet)
# explore_ast(code_snippet)
ast_to_json(code_snippet)