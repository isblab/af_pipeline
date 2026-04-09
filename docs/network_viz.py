# NOTE: portions of this code are taken from github copilot.
import os
import ast
import glob
import json
import getpass
import argparse
import networkx as nx
from collections import defaultdict
from pyvis.network import Network

_user = getpass.getuser()

GITHUB_ORGANIZATION = "isblab"

IGNORE_SUBMODULES = [
    f"/home/{_user}/Projects/af_pipeline/docs",
    f"/home/{_user}/Projects/af_pipeline/tests",
]

NODE_COLORS = {
    "class": "#FFDD00",
    "method": "#00A3FF",
    "function": "#00FF19",
    "unknown": "#D3D3D3",
    "module": "#DA009F",
}

EDGE_COLORS = {
    "defines": "#E3E3E3",
    "calls": "#274c77",
    "inherits": "#ff595e",
    "has_a": "#9d4edd",
    "uses": "#00FF19",
}

def split_with_join(s, delimiter, n):
    """
    Splits the string into two parts at the nth occurrence using split/join.

    Returns a tuple (part_before, part_after).
    """
    parts = s.rsplit(delimiter)
    if len(parts) > n:
        part_before = delimiter.join(parts[:n])
        part_after = delimiter.join(parts[n:])
        return part_before, part_after
    else:
        return s, ""

def create_graph(
    edges: list,
    nodes: dict,
):
    """ Creates a directed graph using NetworkX from the given edges and nodes.

    Args:
        edges (list):
            List of edges where each edge is a dictionary.
            Valid attributes for each edge are:
            - from: The source node name.
            - to: The target node name.
            - type: The type of the edge (e.g., "calls", "defines").
            - color: The color of the edge.
        nodes (dict):
            Dictionary of nodes where each key is a node ID and the value is a dictionary
            containing node attributes. Valid attributes for each node are:
            - name: The name of the node.
            - type: The type of the node (e.g., "class", "method").
            - color: The color of the node.
            - title: The title of the node (used for hover text).

    Returns:
        **G (networkx.classes.digraph.DiGraph)**:
            A directed graph constructed from the provided edges and nodes.
    """

    G = nx.DiGraph()

    for node_id, node_data in nodes.items():

        G.add_node(
        node_id,
        label=node_data["name"],
        type=node_data["type"],
        color=node_data["color"],
        title=node_data.get("title", ""),
        module=node_data.get("module", ""),
        alpha=0.5,
    )

    for edge in edges:
        from_id = next(
            (id for id, data in nodes.items() if data["name"] == edge["from"]), None
        )
        to_id = next(
            (id for id, data in nodes.items() if data["name"] == edge["to"]), None
        )
        if from_id is not None and to_id is not None:
            G.add_edge(
            from_id,
            to_id,
            type=edge["type"],
            label=edge["type"],
            font={"color": edge["color"]},
            color=edge["color"],
            alpha=0.5,
        )

    return G

class CallScopeVisitor(ast.NodeVisitor):
    def __init__(self):
        self.scope_stack = []
        self.orphan_calls = []
        self.method_calls = []

    def visit_ClassDef(self, node):
        self.scope_stack.append(node.name)
        self.generic_visit(node)
        self.scope_stack.pop()

    def visit_FunctionDef(self, node):
        self.scope_stack.append(node.name)
        self.generic_visit(node)
        self.scope_stack.pop()

    def visit_AsyncFunctionDef(self, node):
        self.scope_stack.append(node.name)
        self.generic_visit(node)
        self.scope_stack.pop()

    def visit_Call(self, node):
        # Determine if we are inside a function/class
        if not self.scope_stack:
            call_name = self.get_call_name(node.func)
            if call_name:
                self.orphan_calls.append(call_name)
        else:
            call_name = self.get_call_name(node.func)
            prev_scopes = self.scope_stack[:-1]  # Exclude the current scope
            if call_name:
                self.method_calls.append(f"{'.'.join(prev_scopes)}.{call_name}")

        self.generic_visit(node)

    def get_call_name(self, node):
        """Helper method to extract the function name from different types of nodes."""
        if isinstance(node, ast.Name):
            return node.id  # Simple function call: func()
        elif isinstance(node, ast.Attribute):
            # Method call: obj.method()
            value_name = self.get_call_name(node.value)
            if value_name:
                return f"{value_name}.{node.attr}"
        return None

class CallVisitor(ast.NodeVisitor):
    def __init__(self):
        self.calls = []

    def visit_Call(self, node):
        """Called when a function call (ast.Call) node is encountered."""
        call_name = self.get_call_name(node.func)
        if call_name:
            self.calls.append(call_name)
        # Ensure that child nodes are also visited
        self.generic_visit(node)

    def get_call_name(self, node):
        """Helper method to extract the function name from different types of nodes."""
        if isinstance(node, ast.Name):
            return node.id  # Simple function call: func()
        elif isinstance(node, ast.Attribute):
            # Method call: obj.method()
            value_name = self.get_call_name(node.value)
            if value_name:
                return f"{value_name}.{node.attr}"
        return None

    def get_assigned_name(self, node):
        """Helper method to extract the variable name from an assignment target."""
        if isinstance(node, ast.Name):
            return node.id  # Simple variable assignment: var = ...
        elif isinstance(node, ast.Attribute):
            # Attribute assignment: obj.attr = ...
            value_name = self.get_assigned_name(node.value)
            if value_name:
                return f"{value_name}.{node.attr}"
        return None

class AttributeVisitor(ast.NodeVisitor):
    def __init__(self):
        self.attributes = []

    def visit_Attribute(self, node):
        """Called when an attribute access (ast.Attribute) node is encountered."""
        attr_name = self.get_attribute_name(node)
        if attr_name:
            self.attributes.append(attr_name)
        # Ensure that child nodes are also visited
        self.generic_visit(node)

    def get_attribute_name(self, node):
        """Helper method to extract the full attribute name from nested attributes."""
        if isinstance(node, ast.Name):
            return node.id  # Base case: variable name
        elif isinstance(node, ast.Attribute):
            value_name = self.get_attribute_name(node.value)
            if value_name:
                return f"{value_name}.{node.attr}"
        return None

class LineNoVisitor(ast.NodeVisitor):
    def __init__(self):
        self.func_linenos = {}
        self.class_linenos = {}
        self.class_method_linenos = {}

    def visit_FunctionDef(self, node):
        self.func_linenos[node.name] = node.lineno
        self.generic_visit(node)

    def visit_ClassDef(self, node):
        self.class_linenos[node.name] = node.lineno
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                self.class_method_linenos[f"{node.name}.{item.name}"] = item.lineno
        self.generic_visit(node)

    def visit_AsyncFunctionDef(self, node):
        self.func_linenos[node.name] = node.lineno
        self.generic_visit(node)

class Visitor(ast.NodeVisitor):
    def __init__(self):
        self.functions = []
        self.function_calls = []
        self.classes = []
        self.class_functions = []
        self.class_calls = []
        self.class_inheritances = []
        self.imports = []
        self.calls = []
        self.orphan_calls = []
        self.assignments = []

    def visit_FunctionDef(self, node):
        """Called when a function definition (ast.FunctionDef) node is encountered."""
        self.functions.append(node.name)
        call_visitor = CallVisitor()
        call_visitor.visit(node)
        self.function_calls.append({node.name: call_visitor.calls})
        # Ensure that child nodes are also visited
        self.generic_visit(node)

    def visit_ClassDef(self, node):
        """Called when a class definition (ast.ClassDef) node is encountered."""
        self.classes.append(node.name)
        call_visitor = CallVisitor()
        call_visitor.visit(node)
        self.class_calls.append({node.name: call_visitor.calls})
        self.class_functions.append(
            {node.name: [func.name for func in node.body if isinstance(func, ast.FunctionDef)]}
        )
        self.class_inheritances.append(
            {node.name: [base.id for base in node.bases if isinstance(base, ast.Name)]}
        )
        # Ensure that child nodes are also visited
        self.generic_visit(node)

    def visit_Import(self, node):
        """Called when an import statement (ast.Import) node is encountered."""
        for alias in node.names:
            self.imports.append(alias.name)
        self.generic_visit(node)

    def visit_ImportFrom(self, node):
        """Called when an import from statement (ast.ImportFrom) node is encountered."""
        module_name = node.module if node.module else ""
        for alias in node.names:
            self.imports.append(f"{module_name}.{alias.name}")
        self.generic_visit(node)

    def visit_AsyncFunctionDef(self, node):
        """Called when an async function definition (ast.AsyncFunctionDef) node is encountered."""
        self.functions.append(node.name)
        call_visitor = CallVisitor()
        call_visitor.visit(node)
        self.function_calls.append({node.name: call_visitor.calls})
        # Ensure that child nodes are also visited
        self.generic_visit(node)

    def visit_Call(self, node):
        """Called when a function call (ast.Call) node is encountered."""
        call_visitor = CallVisitor()
        call_visitor.visit(node)
        self.calls.extend(call_visitor.calls)

        call_scope_visitor = CallScopeVisitor()
        call_scope_visitor.visit(node)
        self.orphan_calls.extend(call_scope_visitor.orphan_calls)

        # Ensure that child nodes are also visited
        self.generic_visit(node)

    def visit_Assign(self, node):
        """Called when an assignment (ast.Assign) node is encountered."""
        for target in node.targets:
            assigned_name = CallVisitor().get_assigned_name(target)
            if assigned_name:

                if isinstance(node.value, ast.Call):
                    call_visitor = CallVisitor()
                    call_visitor.visit(node.value)
                    self.assignments.append({assigned_name: call_visitor.calls})

                # elif isinstance(node.value, ast.Attribute):
                #     attribute_visitor = AttributeVisitor()
                #     attribute_visitor.visit(node.value)
                #     self.assignments.append({assigned_name: attribute_visitor.attributes})
                # else:
                #     self.assignments.append({assigned_name: node.value})
        self.generic_visit(node)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--analyse_dirs",
        nargs='*',
        default=[f"/home/{_user}/Projects/af_pipeline"],
        help="Directories to analyze for Python files.",
    )

    parser.add_argument(
        "--module_dirs",
        nargs='*',
        default=[f"/home/{_user}/Projects/af_pipeline/af_pipeline"],
        help="Directories corresponding to the modules being analyzed.",
    )

    parser.add_argument(
        "--module_names",
        nargs='*',
        default=["af_pipeline"],
        help="Names of the modules being analyzed.",
    )

    parser.add_argument(
        "--repository_names",
        nargs="*",
        default=["af_pipeline:af_pipeline"],
        help="Names of the repositories being analyzed in the format 'module_name:'.",
    )

    parser.add_argument(
        "--repository_organization",
        type=str,
        default=GITHUB_ORGANIZATION,
        help="GitHub organization name for constructing repository URLs.",
    )

    parser.add_argument(
        "--only_uses",
        action='store_true',
        default=False,
        help="Only include 'uses' relationships in the graph.",
    )

    parser.add_argument(
        "--output_path",
        type=str,
        default=f"/home/{_user}/Projects/af_pipeline/docs/network/af_pipeline_network.html",
        help="Path to save the output graph html file.",
    )

    parser.add_argument(
        "--show_physics_button",
        action='store_true',
        default=False,
        help="Show physics button in the graph visualization.",
    )

    args = parser.parse_args()

    only_uses = args.only_uses
    analyse_dirs = args.analyse_dirs
    module_dirs = args.module_dirs
    allowed_modules = args.module_names
    repo_names = {
        m_: repo_name
        for m_, repo_name in (repo.split(":") for repo in args.repository_names)
    }
    output_path = os.path.abspath(args.output_path)
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)

    path_prefixes = [os.path.dirname(d) for d in analyse_dirs]

    ############################################################################
    # Recursively look for all the python files in the analyse directories
    all_python_files = [
        (p, f, m) for p, d, m in zip(path_prefixes, analyse_dirs, module_dirs)
        for f in glob.glob(os.path.join(d, "**/*.py"), recursive=True)
    ]
    all_python_files = [(p, f, m) for (p, f, m) in all_python_files]
    all_python_files = [
        (p, f, m) for (p, f, m) in all_python_files
        if not any(f.startswith(ignore) for ignore in IGNORE_SUBMODULES)
    ]
    only_python_file_names = [f for _, f, _ in all_python_files]

    print(f"Found {len(all_python_files)} Python files in the project.")
    ############################################################################

    module_dependencies = defaultdict(set)
    module_calls = defaultdict(set)
    module_orphan_calls = defaultdict(list)
    module_functions = defaultdict(set)
    module_function_calls = defaultdict(dict)
    module_classes = defaultdict(set)
    module_class_calls = defaultdict(dict)
    module_class_functions = defaultdict(dict)
    module_assignments = defaultdict(list)
    module_inheritances = defaultdict(list)

    func_linenos = defaultdict(dict)
    class_linenos = defaultdict(dict)

    # double imports are those modules which have module.module structure
    double_import = {}
    module_file_paths = defaultdict(list)

    for prefix, file_path, module_dir in all_python_files:

        if file_path.startswith(prefix):

            sub_module_name = (
                file_path.replace(prefix, "").replace(".py", "").replace("/", ".").lstrip(".")
            )
            module_file_paths[sub_module_name].append(file_path)

            if os.path.basename(os.path.dirname(module_dir)) == sub_module_name.split(".")[1]:
                double_import[sub_module_name] = True
            else:
                # print(sub_module_name)
                double_import[sub_module_name] = False

    # extend double_import to all hierarchies of modules
    extended_double_import = {}
    for module in double_import:
        parts = module.split(".")
        for i in range(1, len(parts) + 1):
            sub_module = ".".join(parts[:i])
            extended_double_import[sub_module] = double_import[module]

    double_import = extended_double_import

    for module, file_paths in module_file_paths.items():

        for file_path in file_paths:

            with open(file_path, "r") as f:
                content = f.read()

            content_ast = ast.parse(content)

            lineno_visitor = LineNoVisitor()
            lineno_visitor.visit(content_ast)

            for func_name, lineno in lineno_visitor.func_linenos.items():
                full_func_name = f"{module}.{func_name}"
                func_linenos[module][full_func_name] = {
                    "lineno": lineno,
                    "path": file_path,
                }

            for class_name, lineno in lineno_visitor.class_linenos.items():
                full_class_name = f"{module}.{class_name}"
                class_linenos[module][full_class_name] = {
                    "lineno": lineno,
                    "path": file_path,
                }

            for class_method_name, lineno in lineno_visitor.class_method_linenos.items():
                full_class_method_name = f"{module}.{class_method_name}"
                func_linenos[module][full_class_method_name] = {
                    "lineno": lineno,
                    "path": file_path,
                }

            visitor = Visitor()
            visitor.visit(content_ast)

            for imported_ in visitor.imports:
                imported_parent = imported_.rsplit(".", 1)[0]
                if double_import.get(f"{imported_.split(".")[0]}.{imported_parent}", False) is True:
                    module_dependencies[module].add(f"{imported_.split(".")[0]}.{imported_}")
                else:
                    module_dependencies[module].add(imported_)

            for _class in visitor.classes:
                module_classes[module].add(_class)

            for class_func_dict in visitor.class_functions:
                module_class_functions[module].update(class_func_dict)

            for class_call_dict in visitor.class_calls:
                module_class_calls[module].update(class_call_dict)

            for func_call_dict in visitor.function_calls:
                module_function_calls[module].update(func_call_dict)

            for func in visitor.functions:
                _curr_class_funcs = [
                    f for class_func_dict in visitor.class_functions
                    for funcs in class_func_dict.values() for f in funcs
                ]
                if func not in _curr_class_funcs:
                    module_functions[module].add(func)

            for call in visitor.calls:
                module_calls[module].add(call)

            for orphan_call in visitor.orphan_calls:
                module_orphan_calls[module].append(orphan_call)

            for assignment_dict in visitor.assignments:
                module_assignments[module].append(assignment_dict)

            for inheritance_dict in visitor.class_inheritances:
                module_inheritances[module].append(inheritance_dict)

    # this expands the class names to full names
    # format is {module: [class1, class2, ...]}
    module_classes = {
        k: [f"{k}.{cls}" for cls in v]
        for k, v in module_classes.items()
        if k.split(".")[0] in allowed_modules
    }

    # these are class inheritances, which indicate which classes inherit from
    # which other classes.
    # format is {module: [{class_name: [base_class1, base_class2, ...]}, ...]}
    module_inheritances = {
        k: [
            {f"{k}.{class_name}": [f"{base_class}" for base_class in base_classes]}
            for inheritance_dict in v for class_name, base_classes in inheritance_dict.items()
        ] for k, v in module_inheritances.items()
    }

    updated_module_inheritances = defaultdict(list)

    for k, v in module_inheritances.items():
        module_dep_ = module_dependencies.get(k, set())
        for inheritance_dict in v:
            for class_name, base_classes in inheritance_dict.items():
                for base_class in base_classes:
                    if base_class in [x.split(".")[-1] for x in module_dep_]:
                        full_base_class_name = [x for x in module_dep_ if x.split(".")[-1] == base_class][0]
                        if full_base_class_name.split(".")[0] in allowed_modules:
                            updated_module_inheritances[k].append({class_name: [full_base_class_name]})
                    elif base_class in [x.split(".")[-1] for x in module_classes.get(k, set())]:
                        full_base_class_name = [x for x in module_classes.get(k, set()) if x.split(".")[-1] == base_class][0]
                        updated_module_inheritances[k].append({class_name: [full_base_class_name]})

    module_inheritances = updated_module_inheritances

    # this includes all imports, including those from external libraries.
    # We will filter these out later to focus on internal dependencies.
    # this expands the names to full names of the modules
    module_dependencies = {
        k: [dep for dep in v if dep.split(".")[0] in allowed_modules]
        for k, v in module_dependencies.items()
        if k.split(".")[0] in allowed_modules
    }

    # these are orphan functions that are not defined in any internal class.
    # We will filter these later to focus on internal dependencies.
    # this expands the function names to full names
    module_functions = {
        k: [f"{k}.{func}" for func in v]
        for k, v in module_functions.items()
        if k.split(".")[0] in allowed_modules
    }

    # these are class functions (methods), which are defined within a class.
    # format is {module: {class_name: [method1, method2, ...]}}
    # this expands the method names and its parent classes to full names
    module_class_functions = {
        k: {
            f"{k}.{class_name}": [f"{k}.{class_name}.{func}" for func in calls]
            for class_name, calls in v.items()
        } for k, v in module_class_functions.items()
        if k.split(".")[0] in allowed_modules
    }

    flattened_module = (
        list(module_dependencies.keys()) +
        [dep for deps in module_dependencies.values() for dep in deps] +
        list(module_functions.keys()) +
        [func for funcs in module_functions.values() for func in funcs] +
        list(module_classes.keys()) +
        [class_ for classes in module_classes.values() for class_ in classes] +
        list(module_class_functions.keys()) +
        [func for class_funcs in module_class_functions.values() for funcs in class_funcs.values() for func in funcs]
    )

    flattened_module = list(set(flattened_module))

    func_linenos = {
        k: {
            func_name: info for func_name, info in v.items()
            if func_name in flattened_module
        } for k, v in func_linenos.items()
    }

    class_linenos = {
        k: {
            class_name: info for class_name, info in v.items()
            if class_name in flattened_module
        } for k, v in class_linenos.items()
    }

    flat_func_linenos = {
        func_name: info for v in func_linenos.values() for func_name, info in v.items()
    }

    flat_class_linenos = {
        class_name: info for v in class_linenos.values() for class_name, info in v.items()
    }

    # these are calls made within function definitions, which could be calls to
    # methods of the same class or calls to other functions/classes.
    # format is {module: {function_name: [call1, call2, ...]}}
    module_function_calls = {
        k: {
            f"{k}.{func_name}": [
                call for call in calls
                if call in [x.split(".")[-1] for x in module_dependencies.get(k, [])]
                or call in [x.split(".")[-1] for x in module_functions.get(k, [])]
                or call in [x.split(".")[-1] for x in module_classes.get(k, [])]
            ] for func_name, calls in v.items()
            if f"{k}.{func_name}" in flattened_module
        } for k, v in module_function_calls.items()
    }

    updated_module_function_calls = defaultdict(dict)

    for k, v in module_function_calls.items():
        for func_name, calls in v.items():
            for call in calls:

                if call in [x.split(".")[-1] for x in module_dependencies.get(k, [])]:
                    full_call_name = [x for x in module_dependencies.get(k, []) if x.split(".")[-1] == call][0]
                    updated_module_function_calls[k].setdefault(func_name, []).append(full_call_name)

                elif call in [x.split(".")[-1] for x in module_functions.get(k, [])]:
                    full_call_name = [x for x in module_functions.get(k, []) if x.split(".")[-1] == call][0]
                    updated_module_function_calls[k].setdefault(func_name, []).append(full_call_name)

                elif call in [x.split(".")[-1] for x in module_classes.get(k, [])]:
                    full_call_name = [x for x in module_classes.get(k, []) if x.split(".")[-1] == call][0]
                    updated_module_function_calls[k].setdefault(func_name, []).append(full_call_name)

    module_function_calls = updated_module_function_calls

    # these are calls made within class definitions, which could be calls to
    # methods of the same class or calls to other functions/classes.
    # format is {module: {class_name: [call1, call2, ...]}}
    module_class_calls = {
        k: {
            f"{k}.{class_name}": [
                call for call in calls
                if call in [x.split(".")[-1] for x in module_dependencies.get(k, [])]
                or call in [x.split(".")[-1] for x in module_functions.get(k, [])]
                or call in [x.split(".")[-1] for x in module_classes.get(k, [])]
            ] for class_name, calls in v.items()
        } for k, v in module_class_calls.items()
    }

    updated_module_class_calls = defaultdict(dict)

    for k, v in module_class_calls.items():
        for class_name, calls in v.items():
            for call in calls:

                if call in [x.split(".")[-1] for x in module_dependencies.get(k, [])]:
                    full_call_name = [x for x in module_dependencies.get(k, []) if x.split(".")[-1] == call][0]
                    updated_module_class_calls[k].setdefault(class_name, []).append(full_call_name)

                elif call in [x.split(".")[-1] for x in module_functions.get(k, [])]:
                    full_call_name = [x for x in module_functions.get(k, []) if x.split(".")[-1] == call][0]
                    updated_module_class_calls[k].setdefault(class_name, []).append(full_call_name)

                elif call in [x.split(".")[-1] for x in module_classes.get(k, [])]:
                    full_call_name = [x for x in module_classes.get(k, []) if x.split(".")[-1] == call][0]
                    updated_module_class_calls[k].setdefault(class_name, []).append(full_call_name)

    module_class_calls = updated_module_class_calls


    updated_module_assignments = defaultdict(list)
    module_calls = defaultdict(set)
    sequential_full_call_cache = defaultdict(dict)

    for k, v in module_assignments.items():
        for assignment in v:
            for var_name, calls in assignment.items():
                for call in calls:

                    if call in [x.split(".")[-1] for x in module_dependencies.get(k, [])]:
                        full_call_name = [x for x in module_dependencies.get(k, []) if x.split(".")[-1] == call][0]
                        sequential_full_call_cache[k][var_name] = full_call_name
                        updated_module_assignments[k].append({var_name: full_call_name})
                        module_calls[k].add(full_call_name)

                    elif call in [x.split(".")[-1] for x in module_functions.get(k, [])]:
                        full_call_name = [x for x in module_functions.get(k, []) if x.split(".")[-1] == call][0]
                        sequential_full_call_cache[k][var_name] = full_call_name
                        updated_module_assignments[k].append({var_name: full_call_name})
                        module_calls[k].add(full_call_name)

                    elif call in [x.split(".")[-1] for x in module_classes.get(k, [])]:
                        full_call_name = [x for x in module_classes.get(k, []) if x.split(".")[-1] == call][0]
                        sequential_full_call_cache[k][var_name] = full_call_name
                        updated_module_assignments[k].append({var_name: full_call_name})
                        module_calls[k].add(full_call_name)

                    else:
                        call_0 = call.split(".")[0]
                        if call_0 in sequential_full_call_cache[k]:
                            full_call_name = sequential_full_call_cache[k][call_0] + "." + ".".join(call.split(".")[1:])
                            sequential_full_call_cache[k][var_name] = full_call_name
                            updated_module_assignments[k].append({var_name: full_call_name})
                            if full_call_name in flattened_module:
                                module_calls[k].add(full_call_name)

    module_assignments = updated_module_assignments

    # these are internal functions or classes that are called but not defined in
    # any of the internal functions or classes.
    updated_module_orphan_calls = defaultdict(set)

    for k, orphan_calls in module_orphan_calls.items():
        for call in orphan_calls:

            call_0 = call.split(".")[0]

            if call_0 not in sequential_full_call_cache[k]:
                continue
            full_call_name = sequential_full_call_cache[k][call_0] + "." + ".".join(call.split(".")[1:])

            if full_call_name in flattened_module:
                updated_module_orphan_calls[k].add(full_call_name)

    module_orphan_calls = updated_module_orphan_calls

    ############################################################################
    # Creating and visualizing the graph
    ############################################################################

    nodes = {}
    edges = []

    only_modules = (
        set(module_dependencies.keys()) |
        set(module_classes.keys()) |
        set(module_functions.keys()) |
        set(module_class_functions.keys())
    )

    module_hiers = []

    for mod in only_modules:
        _dot_count = mod.count(".")
        temp_mod = mod
        for i in range(_dot_count):
            parent_module = split_with_join(temp_mod, ".", -1)[0]
            if parent_module not in nodes:
                nodes[parent_module] = {
                    "type": "module",
                    "module": parent_module,
                    "name": parent_module,
                    "color": NODE_COLORS["module"],
                    "lineno": None,
                    "file_path": os.path.dirname(module_file_paths.get(parent_module, [""])[0])
                }
            if temp_mod not in nodes:
                nodes[temp_mod] = {
                    "type": "module",
                    "module": temp_mod,
                    "name": temp_mod,
                    "color": NODE_COLORS["module"],
                }
            if temp_mod != parent_module:
                edges.append({
                    "from": parent_module,
                    "to": temp_mod,
                    "type": "contains",
                    "color": EDGE_COLORS["defines"],
                    "font": {"color": EDGE_COLORS["defines"]},
                })
            temp_mod = parent_module

    if only_uses is True:

        for k, module_assigns in module_assignments.items():
            for assign_dict in module_assigns:
                for var_name, call in assign_dict.items():
                    if call in [class_ for classes in module_classes.values() for class_ in classes]:
                        call_type = "class"
                        color_ = NODE_COLORS["class"]
                        lineno = class_linenos.get(k, {}).get(call, {}).get("lineno", None)
                        file_path = class_linenos.get(k, {}).get(call, {}).get("path", "")
                    elif (
                        call in [func_ for funcs in module_functions.values() for func_ in funcs] or
                        call in [func_ for class_funcs in module_class_functions.values() for funcs in class_funcs.values() for func_ in funcs]
                    ):
                        call_type = "function"
                        color_ = NODE_COLORS["function"]
                        lineno = func_linenos.get(k, {}).get(call, {}).get("lineno", None)
                        file_path = func_linenos.get(k, {}).get(call, {}).get("path", "")

                    if call not in nodes:
                        nodes[call] = {
                            "type": call_type,
                            "module": k,
                            "name": call,
                            "color": color_,
                            "lineno": lineno,
                            "file_path": file_path,
                        }
                    if k not in nodes:
                        nodes[k] = {
                            "type": "module",
                            "module": k,
                            "name": k,
                            "color": NODE_COLORS["module"],
                            "lineno": None,
                            "file_path": os.path.dirname(module_file_paths.get(k, [None])[0]) if module_file_paths.get(k, [None])[0] else None,
                        }
                    edges.append({
                        "from": k,
                        "to": call,
                        "type": "uses",
                        "color": EDGE_COLORS["uses"],
                        "font": {"color": EDGE_COLORS["uses"]},
                    })

    elif only_uses is False:

        for k, class_func_dict in module_class_functions.items():
            for class_name, funcs in class_func_dict.items():
                nodes[class_name] = {
                    "type": "class",
                    "module": k,
                    "name": class_name,
                    "color": NODE_COLORS["class"],
                    "lineno": class_linenos.get(k, {}).get(class_name, None).get("lineno", None),
                    "file_path": class_linenos.get(k, {}).get(class_name, None).get("path", None),
                }
                nodes[k] = {
                    "type": "module",
                    "module": k,
                    "name": k,
                    "color": NODE_COLORS["module"],
                    "lineno": None,
                    "file_path": os.path.dirname(module_file_paths.get(k, [None])[0]) if module_file_paths.get(k, [None])[0] else None,
                }
                edges.append({
                    "from": k,
                    "to": class_name,
                    "type": "defines",
                    "color": EDGE_COLORS["defines"],
                    "font": {"color": EDGE_COLORS["defines"]},
                })
                for func in funcs:
                    nodes[func] = {
                        "type": "function",
                        "module": k,
                        "name": func,
                        "color": NODE_COLORS["function"],
                        "lineno": func_linenos.get(k, {}).get(func, None).get("lineno", None),
                        "file_path": func_linenos.get(k, {}).get(func, None).get("path", None),
                    }
                    edges.append({
                        "from": class_name,
                        "to": func,
                        "type": "defines",
                        "color": EDGE_COLORS["defines"],
                        "font": {"color": EDGE_COLORS["defines"]},
                    })

        for k, orphan_calls in module_orphan_calls.items():
            for call in orphan_calls:
                if call in [class_ for classes in module_classes.values() for class_ in classes]:
                    call_type = "class"
                    color_ = NODE_COLORS["class"]
                    lineno = class_linenos.get(k, {}).get(call, {}).get("lineno", None)
                    file_path = class_linenos.get(k, {}).get(call, {}).get("path", None)
                elif (
                    call in [func_ for funcs in module_functions.values() for func_ in funcs] or
                    call in [func_ for class_funcs in module_class_functions.values() for funcs in class_funcs.values() for func_ in funcs]
                ):
                    call_type = "function"
                    color_ = NODE_COLORS["function"]
                    lineno = flat_func_linenos.get(call, {}).get("lineno", None)
                    file_path = flat_func_linenos.get(call, {}).get("path", None)
                else:
                    call_type = "unknown"
                    color_ = NODE_COLORS["unknown"]
                    lineno = None
                    file_path = os.path.dirname(module_file_paths.get(k, [""])[0])
                if call not in nodes:
                    nodes[call] = {
                        "type": call_type,
                        "module": k,
                        "name": call,
                        "color": color_,
                        "lineno": lineno,
                        "file_path": file_path,
                    }
                    nodes[k] = {
                        "type": "module",
                        "module": k,
                        "name": k,
                        "color": NODE_COLORS["module"],
                        "lineno": None,
                        "file_path": os.path.dirname(module_file_paths.get(k, [None])[0]) if module_file_paths.get(k, [None])[0] else None,
                    }
                    edges.append({
                        "from": k,
                        "to": call,
                        "type": "calls",
                        "color": EDGE_COLORS["calls"],
                        "font": {"color": EDGE_COLORS["calls"]},
                    })

        for k, func_call_dict in module_function_calls.items():
            for func_name, calls in func_call_dict.items():
                for call in calls:
                    if call in [class_ for classes in module_classes.values() for class_ in classes]:
                        call_type = "class"
                        color_ = NODE_COLORS["class"]
                        lineno = flat_class_linenos.get(call, {}).get("lineno", None)
                        file_path = flat_class_linenos.get(call, {}).get("path", None)
                    elif (
                        call in [func_ for funcs in module_functions.values() for func_ in funcs] or
                        call in [func_ for class_funcs in module_class_functions.values() for funcs in class_funcs.values() for func_ in funcs]
                    ):
                        call_type = "function"
                        color_ = NODE_COLORS["function"]
                        lineno = flat_func_linenos.get(call, {}).get("lineno", None)
                        file_path = flat_func_linenos.get(call, {}).get("path", None)
                    else:
                        call_type = "unknown"
                        color_ = NODE_COLORS["unknown"]
                        lineno = None
                        file_path = os.path.dirname(module_file_paths.get(k, [None])[0]) if module_file_paths.get(k, [None])[0] else None
                    if call not in nodes:
                        nodes[call] = {
                            "type": call_type,
                            "module": k,
                            "name": call,
                            "color": color_,
                            "lineno": lineno,
                            "file_path": file_path,
                        }
                    edges.append({
                        "from": func_name,
                        "to": call,
                        "type": "calls",
                        "color": EDGE_COLORS["calls"],
                        "font": {"color": EDGE_COLORS["calls"]},
                    })

        for k, class_call_dict in module_class_calls.items():
            for class_name, calls in class_call_dict.items():
                for call in calls:
                    if call in [class_ for classes in module_classes.values() for class_ in classes]:
                        call_type = "class"
                        color_ = NODE_COLORS["class"]
                        relation = "has_a"
                        lineno = flat_class_linenos.get(call, {}).get("lineno", None)
                        file_path = flat_class_linenos.get(call, {}).get("path", None)
                    elif (
                        call in [func_ for funcs in module_functions.values() for func_ in funcs] or
                        call in [func_ for class_funcs in module_class_functions.values() for funcs in class_funcs.values() for func_ in funcs]
                    ):
                        call_type = "function"
                        color_ = NODE_COLORS["function"]
                        relation = "calls"
                        lineno = flat_func_linenos.get(call, {}).get("lineno", None)
                        file_path = flat_func_linenos.get(call, {}).get("path", None)
                    else:
                        call_type = "unknown"
                        color_ = NODE_COLORS["unknown"]
                        relation = "calls"
                        lineno = None
                        file_path = os.path.dirname(module_file_paths.get(k, [None])[0]) if module_file_paths.get(k, [None])[0] else None
                    if call not in nodes:
                        nodes[call] = {
                            "type": call_type,
                            "module": k,
                            "name": call,
                            "color": color_,
                            "lineno": lineno,
                            "file_path": file_path,
                        }
                    if class_name not in nodes:
                        nodes[class_name] = {
                            "type": "class",
                            "module": k,
                            "name": class_name,
                            "color": NODE_COLORS["class"],
                            "lineno": flat_class_linenos.get(class_name, {}).get("lineno", None),
                            "file_path": flat_class_linenos.get(class_name, {}).get("path", None),
                        }
                    edges.append({
                        "from": class_name,
                        "to": call,
                        "type": relation,
                        "color": EDGE_COLORS[relation],
                        "font": {"color": EDGE_COLORS[relation]},
                    })

        for k, module_inheritance_list in module_inheritances.items():
            for inheritance_dict in module_inheritance_list:
                for class_name, base_classes in inheritance_dict.items():
                    for base_class in base_classes:
                        if base_class not in nodes:
                            nodes[base_class] = {
                                "type": "class",
                                "module": k,
                                "name": base_class,
                                "color": NODE_COLORS["class"],
                                "lineno": flat_class_linenos.get(base_class, {}).get("lineno", None),
                                "file_path": flat_class_linenos.get(base_class, {}).get("path", None),
                            }
                        if class_name not in nodes:
                            nodes[class_name] = {
                                "type": "class",
                                "module": k,
                                "name": class_name,
                                "color": NODE_COLORS["class"],
                                "lineno": flat_class_linenos.get(class_name, {}).get("lineno", None),
                                "file_path": flat_class_linenos.get(class_name, {}).get("path", None),
                            }
                        edges.append({
                            "from": class_name,
                            "to": base_class,
                            "type": "inherits",
                            "color": EDGE_COLORS["inherits"],
                            "font": {"color": EDGE_COLORS["inherits"]},
                        })

        for k, funcs_list in module_functions.items():
            for func in funcs_list:
                if func not in nodes:
                    nodes[func] = {
                        "type": "function",
                        "module": k,
                        "name": func,
                        "color": NODE_COLORS["function"],
                        "lineno": func_linenos.get(k, {}).get(func, {}).get("lineno", None),
                        "file_path": func_linenos.get(k, {}).get(func, {}).get("path", None),
                    }
                if k not in nodes:
                    nodes[k] = {
                        "type": "module",
                        "module": k,
                        "name": k,
                        "color": NODE_COLORS["module"],
                        "lineno": None,
                        "file_path": os.path.dirname(module_file_paths.get(k, [None])[0]) if module_file_paths.get(k, [None])[0] else None,
                    }
                edges.append({
                    "from": k,
                    "to": func,
                    "type": "defines",
                    "color": EDGE_COLORS["defines"],
                    "font": {"color": EDGE_COLORS["defines"]},
                })

    G = create_graph(edges, nodes)

    net = Network(
        height="90vh",
        width="98vw",
        directed=True,
        notebook=False,
        layout=False,
        select_menu=True,
        filter_menu=True,
        cdn_resources="in_line",
        # font_color='#10000000',
    )

    # edit node label to show only the basename
    for node_id, data in G.nodes(data=True):

        link_to_add = ""
        node_data = nodes[node_id]
        lineno = node_data.get("lineno", 1)

        if data["type"] not in ["class", "function", "module"]:
            continue

        module_name = data["module"]
        repo_name = repo_names.get(module_name.split(".")[0], None)
        if double_import.get(module_name, False) is True:
            module_base = module_name.split(".")[0]
        else:
            module_base = ""
        if repo_name is None:
            continue

        if data["type"] == "module":

            script_path = module_name.replace(".", "/") + ".py"
            found_python_file = False

            for pref in path_prefixes:
                script_path_ = os.path.join(pref, script_path)
                if script_path_ in only_python_file_names:
                    script_path = script_path_
                    found_python_file = True
                    break

            if found_python_file is False:
                script_path =  module_name.replace(".", "/")

            for idx, module_dir in enumerate(module_dirs):
                script_path = os.path.join(path_prefixes[idx], script_path)
                script_path = script_path.replace(module_dir, "")

            if script_path in analyse_dirs:
                link_to_add = f"https://github.com/{GITHUB_ORGANIZATION}/{repo_name}"
            else:
                for analyse_dir in analyse_dirs:
                    script_path = script_path.replace(analyse_dir, "")
                link_to_add = f"https://github.com/{GITHUB_ORGANIZATION}/{repo_name}/tree/main/{repo_name}/{script_path}"
                if module_base != "":
                    link_to_add = f"https://github.com/{GITHUB_ORGANIZATION}/{repo_name}/tree/main/{module_base}/{script_path}"
                else:
                    # for analyse_dir in analyse_dirs:
                    #     script_path = script_path.replace(analyse_dir, "")
                    link_to_add = f"https://github.com/{GITHUB_ORGANIZATION}/{repo_name}/tree/main/{script_path}"
                    # print(script_path)
                    # print(link_to_add)
        else:
            script_path = node_data.get("file_path", "")
            for analyse_dir in analyse_dirs:
                script_path = script_path.replace(analyse_dir + os.sep, "")
            link_to_add = f"https://github.com/{GITHUB_ORGANIZATION}/{repo_name}/tree/main/{script_path}#L{lineno}"

        data["title"] = (
            f"<a href='{link_to_add}' target='_blank'>{data["label"]}</a>"
        )
        data["label"] = data["label"].split(".")[-1]

    # G.nodes(data=True)
    net.from_nx(G)

    options = {
      "physics": {
        "stabilization": {
          "enabled": True,
          "iterations": 500,
          "updateInterval": 100,
          "fit": True,
        },
        "forceAtlas2Based": {
          "theta": 0.6,
          "gravitationalConstant": -99,
          "springLength": 180,
          "springConstant": 0.15,
          "centralGravity": 0.005,
          "avoidOverlap": 0.8,
        },
        "minVelocity": 0.6,
        "maxVelocity": 100,
        "solver": "forceAtlas2Based"
      },
      "wind": {
        "x": 8.5,
        "y": 0.0,
      },
      "edges": {
        "smooth": {
            "type": "cubicBezier",
            "forceDirection": "vertical",
            "roundedness": 0.6,
        },
      },
      "interaction": {
        "hover": True,
        "multiselect": True,
        "navigationButtons": True,
        "tooltipDelay": 150,
      },
      "layout": {
        "hierarchical": {
        "enabled": False,
        "direction": "UD",
        "sortMethod": "directed",
        "shakeTowards": "leaves",
        "nodeSpacing": 150,
        "levelSeparation": 250,
        "treeSpacing": 80,
        }
      },
    }

    if args.show_physics_button:
        options.update({
            "configure": {
                "enabled": True,
                "filter": "physics",
                "showButton": args.show_physics_button,
            }
        })

    net.set_options(json.dumps(options))
    net.set_template_dir(
        f"/home/{_user}/Projects/af_pipeline/docs/template",
        template_file="template_custom.html"
    )

    # custom javascript to open links in new tab on double click to nodes
    custom_js = """
    <script type="text/javascript">
        network.on("doubleClick", function (params) {
            if (params.nodes.length > 0) {
                var nodeId = params.nodes[0];
                var node = nodes.get(nodeId);
                if (node && node.title) {
                    // parse the HTML in the title to extract any anchor href
                    var wrapper = document.createElement('div');
                    wrapper.innerHTML = node.title;
                    var a = wrapper.querySelector('a');
                    if (a && a.getAttribute('href')) {
                        window.open(a.getAttribute('href'), '_blank');
                    }
                }
            }
        });
    </script>
    """

    net_html = net.generate_html()
    net_html = net_html.replace("</body>", f"{custom_js}</body>")

    os.makedirs(output_dir, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(net_html)
    print(f"Network visualization saved to {output_path}")