import ast
import os
from pprint import pprint
from collections import defaultdict
import networkx as nx
from pyvis.network import Network
import json


def create_graph(edges, nodes):
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
        label=node_data["name"].replace("af_pipeline.", ""),
        type=node_data["type"],
        color=node_data["color"],
        title=node_data.get("title", ""),
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
            color=edge["color"],
        )

    print(type(G))

    return G

def extract_all(script_path, project_dir, module_name=None):
    """
    Extracts class and method names from a Python script.

    Args:
        script_path (str): The path to the Python script.

    Returns:
        **tuple**:
        A tuple containing two lists:
            - class_names (list): Names of classes found in the script.
            - method_names (list): Names of methods found within classes.
    """
    all_classes = []
    class_dict = {}
    all_methods = []
    orphan_funcs = []
    import_dict = defaultdict(dict)
    all_imports = []
    class_inheritance = {}

    curr_module = os.path.splitext(script_path.replace(project_dir, ""))[0].replace(os.path.sep, ".")

    if module_name is None:
        module_name = project_dir.split(os.path.sep)[-1]

    curr_module = f"{module_name}{curr_module}"
    # print(curr_module)

    with open(script_path, "r") as file:
        tree = ast.parse(file.read())

    for node in ast.walk(tree):

        if isinstance(node, ast.ClassDef): # class definition
            all_classes.append(f"{curr_module}.{node.name}")
            class_dict[f"{curr_module}.{node.name}"] = []
            class_inheritance[f"{curr_module}.{node.name}"] = []
            # print(ast.dump(node, indent=4))
            # exit()
            for item in node.body:
                if isinstance(item, ast.FunctionDef): # method definition
                    all_methods.append(f"{curr_module}.{node.name}.{item.name}")
                    class_dict[f"{curr_module}.{node.name}"].append(
                        f"{curr_module}.{node.name}.{item.name}"
                    )
            for base in node.bases: # inheritance
                if isinstance(base, ast.Name):
                    class_inheritance[f"{curr_module}.{node.name}"].append(base.id)
                # elif isinstance(base, ast.Attribute):
                #     parts = []
                #     curr = base
                #     while isinstance(curr, ast.Attribute):
                #         parts.append(curr.attr)
                #         curr = curr.value
                #     if isinstance(curr, ast.Name):
                #         parts.append(curr.id)
                #     parts.reverse()
                #     class_inheritance[f"{curr_module}.{node.name}"].append(".".join(parts))

        elif isinstance(node, ast.FunctionDef): # function definition (includes methods)
            if node.name.startswith("__") is False:
                orphan_funcs.append(f"{curr_module}.{node.name}")

        elif isinstance(node, ast.Import): # imported modules
            for alias in node.names:
                if curr_module.split(".")[0] in alias.name:
                    import_dict[curr_module][alias.name] = [node.lineno, node.end_lineno]

        elif isinstance(node, ast.ImportFrom): # imported modules
            for alias in node.names:
                if curr_module.split(".")[0] in node.module:
                    import_dict[curr_module][f"{node.module}.{alias.name}"] = alias.lineno
                    all_imports.append(f"{node.module}.{alias.name}")

    orphan_funcs = list(set([func for func in orphan_funcs if func not in all_methods]))
    class_inheritance = {k: v for k, v in class_inheritance.items() if v}

    for k, v in class_inheritance.items():
        new_v = []
        for class_name in v:
            candidates = [c for c in all_classes if c.endswith(f".{class_name}")]
            imp_candidates = [c for c in all_imports if c.endswith(f".{class_name}")]
            if candidates:
                # prefer classes defined in same module
                pref = [c for c in candidates if c.startswith(curr_module + ".")]
                selected = pref[0] if pref else candidates[0]
                new_v.append(selected)
            else:
                if imp_candidates:
                    new_v.append(imp_candidates[0])

        class_inheritance[k] = new_v

    return all_classes, class_dict, all_methods, orphan_funcs, import_dict, all_imports, class_inheritance


class DependencyAnalyzer(ast.NodeVisitor):
    def __init__(self, module_dir=None, script_path=None):
        self.module_dir = module_dir
        self.script_path = script_path
        self.module_prefix = self._compute_module_prefix(module_dir, script_path)

        # base package name (used to decide whether an import is "internal")
        self.base_pkg = os.path.basename(module_dir.rstrip(os.sep)) if module_dir else ""

        self.dependencies = {}
        self.current_scope = []  # stack of full names: module.Class, module.Class.method, module.func
        self.var_stack = [dict()]  # stack of variable->full class name mappings for scopes
        self.instance_attrs = defaultdict(dict)  # full_class_name -> {attr_name: full_class_name}
        self.class_names = set()  # set of full class names
        self.simple_name_map = defaultdict(list)  # simple name -> list of full names (classes and functions)

        # import tracking:
        # imported_modules: alias -> full module path (only for internal modules in module_dir)
        # imported_attrs: alias -> full attribute path (e.g. ClassName -> package.module.ClassName) for from-imports (only internal)
        self.imported_modules = {}
        self.imported_attrs = {}

    def _compute_module_prefix(self, module_dir, script_path):
        try:
            rel = os.path.relpath(script_path, module_dir)
        except Exception:
            rel = os.path.basename(script_path)
        parts = rel.split(os.sep)
        # remove __init__.py
        if parts[-1] == "__init__.py":
            parts = parts[:-1]
        else:
            parts[-1] = os.path.splitext(parts[-1])[0]
        base_pkg = os.path.basename(module_dir.rstrip(os.sep))
        if parts and parts[0].startswith(".."):
            # script not inside module_dir; fallback to filename
            return parts[-1]
        if parts:
            return ".".join([base_pkg] + parts) if base_pkg else ".".join(parts)
        return base_pkg or ""

    def _is_internal(self, full_name):
        # Decide whether a module/class path belongs to the provided module_dir package
        if not full_name or not self.base_pkg:
            return False
        return full_name == self.base_pkg or full_name.startswith(self.base_pkg + ".")

    def _is_allowed_target(self, callee):
        """
        Return True if callee should be tracked according to:
         - present in this script (module_prefix)
         - or imported from some internal module (recorded in imported_modules/imported_attrs)
        """
        if not callee:
            return False

        # present in this script
        if self.module_prefix and callee.startswith(self.module_prefix + "."):
            return True

        # imported module prefixes
        for mod in self.imported_modules.values():
            if callee == mod or callee.startswith(mod + "."):
                return True

        # imported attributes (from-imports)
        for attr_full in self.imported_attrs.values():
            if callee == attr_full or callee.startswith(attr_full + "."):
                return True

        return False

    def _push_scope(self):
        self.var_stack.append({})

    def _pop_scope(self):
        self.var_stack.pop()

    def _set_var(self, name, class_full_name):
        # store resolved full class name for variable in current local scope
        if class_full_name:
            self.var_stack[-1][name] = class_full_name

    def _resolve_var(self, name):
        for scope in reversed(self.var_stack):
            if name in scope:
                return scope[name]
        return None

    def _current_class(self):
        # Return the nearest enclosing class full name or None
        for name in reversed(self.current_scope):
            # class full names are registered in class_names
            if name in self.class_names:
                return name
            # also handle nested names like module.Class.method (method entry in current_scope)
            if '.' in name:
                candidate = name.rsplit('.', 1)[0]
                if candidate in self.class_names:
                    return candidate
        return None

    def _register_def(self, kind, simple_name, full_name):
        # kind: 'class' or 'func'
        self.simple_name_map[simple_name].append(full_name)
        if full_name not in self.dependencies:
            self.dependencies[full_name] = set()

    def visit_Import(self, node):
        # import pkg.mod as alias
        for alias in node.names:
            full_mod = alias.name  # e.g. "af_pipeline.submodule" or "requests"
            asname = alias.asname or full_mod.split('.')[0]
            # only keep mapping for internal modules (inside module_dir package)
            if self._is_internal(full_mod):
                self.imported_modules[asname] = full_mod
        self.generic_visit(node)

    def visit_ImportFrom(self, node):
        # from pkg.module import Name as alias
        module = node.module or ""
        level = getattr(node, "level", 0) or 0

        # resolve relative imports against module_prefix more robustly
        if level:
            curr_parts = self.module_prefix.split('.') if self.module_prefix else []
            # remove `level` items from the end of current prefix (approximation)
            if level <= len(curr_parts):
                base_parts = curr_parts[:-level]
            else:
                base_parts = []
            module_parts = module.split('.') if module else []
            full_module = ".".join([p for p in (base_parts + module_parts) if p])
        else:
            full_module = module

        for alias in node.names:
            # skip "from ... import *"
            if alias.name == "*":
                continue

            asname = alias.asname or alias.name

            # if we resolved a full module and it's internal, register imported attr
            if full_module and self._is_internal(full_module):
                full_name = f"{full_module}.{alias.name}"
                # register as a known simple name mapping to an internal full path
                self.simple_name_map[asname].append(full_name)
                self.imported_attrs[asname] = full_name
            else:
                # handle cases like "from . import Name" where full_module was empty
                if level:
                    curr_parts = self.module_prefix.split('.') if self.module_prefix else []
                    if level <= len(curr_parts):
                        base = ".".join(curr_parts[:-level])
                    else:
                        base = ""
                    if base and self._is_internal(base):
                        full_name = f"{base}.{alias.name}"
                        self.simple_name_map[asname].append(full_name)
                        self.imported_attrs[asname] = full_name
                    # else: external or couldn't resolve; ignore for internal mapping

        self.generic_visit(node)

    def visit_ClassDef(self, node):
        # compute full class name based on current scope
        simple_name = node.name
        if self.current_scope:
            parent = self.current_scope[-1]
            # if parent is a class, nest under it
            if parent in self.class_names:
                full_name = f"{parent}.{simple_name}"
            else:
                # parent might be a function or module-level name; treat as top-level in module
                full_name = f"{self.module_prefix}.{simple_name}" if self.module_prefix else simple_name
        else:
            full_name = f"{self.module_prefix}.{simple_name}" if self.module_prefix else simple_name

        self.class_names.add(full_name)
        self._register_def('class', simple_name, full_name)

        # prepare scope as full class name
        self.current_scope.append(full_name)
        self._push_scope()  # class-level scope
        self.generic_visit(node)
        self._pop_scope()
        self.current_scope.pop()

    def visit_FunctionDef(self, node):
        simple_name = node.name
        if self.current_scope:
            parent = self.current_scope[-1]
            # If inside a class, make it a method
            if parent in self.class_names:
                full_name = f"{parent}.{simple_name}"
            else:
                # nested function: qualify with parent full name
                full_name = f"{parent}.{simple_name}"
        else:
            # module-level function
            full_name = f"{self.module_prefix}.{simple_name}" if self.module_prefix else simple_name

        self._register_def('func', simple_name, full_name)

        self.current_scope.append(full_name)
        self._push_scope()  # function/method local scope
        # also treat arguments as local variables (not typed)
        for arg in getattr(node.args, "args", []):
            if isinstance(arg, ast.arg):
                self._set_var(arg.arg, None)
        self.generic_visit(node)
        self._pop_scope()
        self.current_scope.pop()

    def visit_Assign(self, node):
        # Track simple patterns: var = ClassName() and self.attr = ClassName()
        if isinstance(node.value, ast.Call):
            # func can be Name or Attribute (e.g., module.Class())
            target_class_full = None
            if isinstance(node.value.func, ast.Name):
                class_simple = node.value.func.id
                # try to resolve to a full class name using simple_name_map or imported attributes
                candidates = self.simple_name_map.get(class_simple, [])
                if len(candidates) == 1:
                    target_class_full = candidates[0]
                elif len(candidates) > 1:
                    # prefer classes defined in same module
                    pref = [c for c in candidates if c.startswith(self.module_prefix + ".")]
                    target_class_full = pref[0] if pref else candidates[0]
                else:
                    # check imported attrs (from-imports)
                    imported = self.imported_attrs.get(class_simple)
                    if imported:
                        target_class_full = imported
                    else:
                        # not found among definitions: qualify with module_prefix if plausible
                        if self.module_prefix:
                            target_class_full = f"{self.module_prefix}.{class_simple}"
                        else:
                            target_class_full = class_simple

            elif isinstance(node.value.func, ast.Attribute):
                # e.g., pkg.ClassName() -> try to reconstruct dotted name; then try to expand imported module alias
                parts = []
                f = node.value.func
                while isinstance(f, ast.Attribute):
                    parts.append(f.attr)
                    f = f.value
                if isinstance(f, ast.Name):
                    parts.append(f.id)
                parts.reverse()
                candidate = ".".join(parts)
                # if first token corresponds to an imported internal module alias, expand it
                tokens = candidate.split('.')
                if tokens and tokens[0] in self.imported_modules:
                    tokens[0] = self.imported_modules[tokens[0]]
                    candidate = ".".join(tokens)
                target_class_full = candidate

            for target in node.targets:
                if isinstance(target, ast.Name):
                    # variable assignment
                    self._set_var(target.id, target_class_full)
                elif isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name) and target.value.id == 'self':
                    parent_class = self._current_class()
                    if parent_class:
                        self.instance_attrs[parent_class][target.attr] = target_class_full
        self.generic_visit(node)

    def visit_Call(self, node):
        if not self.current_scope:
            self.generic_visit(node)
            return

        caller = self.current_scope[-1]
        callee = None

        # Direct calls like foo()
        if isinstance(node.func, ast.Name):
            name = node.func.id
            # Prefer local definitions/imported attrs
            candidates = self.simple_name_map.get(name, [])
            if len(candidates) == 1:
                callee = candidates[0]
            elif len(candidates) > 1:
                pref = [c for c in candidates if c.startswith(self.module_prefix + ".")]
                callee = pref[0] if pref else candidates[0]
            else:
                # check imported from-X (explicit attr imports)
                imported = self.imported_attrs.get(name)
                if imported:
                    callee = imported
                else:
                    # not found; qualify with module prefix if present (likely a module-level function)
                    callee = f"{self.module_prefix}.{name}" if self.module_prefix else name

        # Attribute calls like obj.method(), self.method(), self.attr.method()
        elif isinstance(node.func, ast.Attribute):
            value = node.func.value

            # obj.method() where obj is a variable: try to resolve its type
            if isinstance(value, ast.Name):
                var_name = value.id
                if var_name == 'self':
                    parent_class = self._current_class()
                    if parent_class:
                        callee = f"{parent_class}.{node.func.attr}"
                    else:
                        callee = f"{self.module_prefix}.{node.func.attr}" if self.module_prefix else node.func.attr
                else:
                    # first try variable instance mapping (var -> class instance)
                    type_name = self._resolve_var(var_name)
                    if type_name:
                        callee = f"{type_name}.{node.func.attr}"
                    else:
                        # maybe var_name is actually a class name (static method call) imported or defined
                        candidates = self.simple_name_map.get(var_name, [])
                        class_candidate = None
                        if candidates:
                            # prefer classes from this module or entries that are known classes
                            pref = [c for c in candidates if c in self.class_names or c.startswith(self.module_prefix + ".")]
                            if pref:
                                # prefer actual class names
                                for c in pref:
                                    if c in self.class_names:
                                        class_candidate = c
                                        break
                                if not class_candidate:
                                    class_candidate = pref[0]
                            else:
                                class_candidate = candidates[0]
                        if class_candidate:
                            callee = f"{class_candidate}.{node.func.attr}"
                        else:
                            # check imported_attrs (from-imported classes)
                            imported = self.imported_attrs.get(var_name)
                            if imported:
                                callee = f"{imported}.{node.func.attr}"
                            else:
                                # fallback: maybe module-level function or external
                                callee = f"{self.module_prefix}.{node.func.attr}" if self.module_prefix else node.func.attr

            # self.attr.method() -> resolve self.attr from instance_attrs if possible
            elif isinstance(value, ast.Attribute) and isinstance(value.value, ast.Name) and value.value.id == 'self':
                attr_name = value.attr
                parent_class = self._current_class()
                if parent_class:
                    mapped = self.instance_attrs.get(parent_class, {}).get(attr_name)
                    if mapped:
                        callee = f"{mapped}.{node.func.attr}"
                    else:
                        callee = f"{self.module_prefix}.{node.func.attr}" if self.module_prefix else node.func.attr
                else:
                    callee = f"{self.module_prefix}.{node.func.attr}" if self.module_prefix else node.func.attr

            else:
                # fallback: try to construct dotted name for call like pkg.obj.method() or Class.method()
                parts = []
                f = node.func
                while isinstance(f, ast.Attribute):
                    parts.append(f.attr)
                    f = f.value
                if isinstance(f, ast.Name):
                    parts.append(f.id)
                parts.reverse()
                candidate = ".".join(parts)
                tokens = candidate.split('.')
                # if the first token is an imported internal module alias, expand it
                if tokens and tokens[0] in self.imported_modules:
                    tokens[0] = self.imported_modules[tokens[0]]
                    candidate = ".".join(tokens)
                # if the first token is a from-imported attr (class), expand it
                if tokens and tokens[0] in self.imported_attrs:
                    base = self.imported_attrs[tokens[0]]
                    rest = tokens[1:]
                    candidate = ".".join([base] + rest) if rest else base
                callee = candidate

        # Add dependency only if resolved and either:
        #  - the callee is defined in this script (module_prefix)
        #  - OR the callee refers to something explicitly imported from an internal module
        if callee:
            if self._is_allowed_target(callee):
                if caller not in self.dependencies:
                    self.dependencies[caller] = set()
                if callee != caller:
                    self.dependencies[caller].add(callee)

        self.generic_visit(node)

def get_dependency_graph(script_path, module_dir):

    with open(script_path, "r") as file:
        tree = ast.parse(file.read())
    analyzer = DependencyAnalyzer(module_dir=module_dir, script_path=script_path)
    analyzer.visit(tree)
    # convert sets to lists for nicer printing if needed
    return {k: set(v) for k, v in analyzer.dependencies.items()}

def get_type(name, classes, methods, orphan_funcs, all_imports):
    if name in classes:
        return "class"
    elif name in methods:
        return "method"
    elif name in orphan_funcs:
        return "function"
    # elif name in all_imports:
    #     return "import"
    else:
        return "unknown"

NODE_COLORS = {
    "class": "#FFDD00",
    "method": "#00A3FF",
    "function": "#00FF19",
    # "import": "#FF00F6",
    "unknown": "#D3D3D3",
}

EDGE_COLORS = {
    "defines": "#D3D3D3",
    "calls": "#274c77",
    "inherits": "#ff595e",
    "has_a": "#09AAc8",
}

if __name__ == "__main__":

    module_dir = "/home/omkar/Projects/af_pipeline/af_pipeline"
    # script_path = "/home/omkar/Projects/af_pipeline/af_pipeline/parser1/structure_parser1.py"
    # script_path = "/home/omkar/Projects/af_pipeline/af_pipeline/parser1/data_parser.py"
    # script_path = "/home/omkar/Projects/af_pipeline/af_pipeline/af_input/alphafold3.py"
    # script_path = "/home/omkar/Projects/af_pipeline/af_pipeline/rank_predictions/rank_af.py"
    # script_path = "/home/omkar/Projects/af_pipeline/af_pipeline/rigid_bodies/rigid_bodies.py"

    # classes, class_dict, methods, orphan_funcs, imports, all_imports, class_inheritance = extract_all(script_path, module_dir)

    # dependency_graph = get_dependency_graph(script_path, module_dir)

    # filtered_graph = {}

    # for k, v in dependency_graph.items():
    #     filtered_graph[k] = set()
    #     for dep in v:
    #         if dep in classes or dep in methods or dep in orphan_funcs or dep in all_imports:
    #             filtered_graph[k].add(dep)

    # filtered_graph = {k: v for k, v in filtered_graph.items() if v}

    # pprint(filtered_graph)
    # pprint(class_dict)
    # pprint(class_inheritance)
    # pprint(all_imports)

    module_dir = "/home/omkar/Projects/af_pipeline/af_pipeline"
    # module_dir = "/home/omkar/Projects/cardiac_desmosome"
    # module_name = "af_pipeline"

    ignore_dirs = [
        # "af_input",
        # "constants",
        # "pae_to_domains",
        # "parser1",
        # "rank_predictions",
        # "rigid_bodies",
        "tests",
        # "tools",
        # "utils",
    ]

    ignore_files = [
        # "_initialize.py",
        # "alphafold3.py",
    ]

    script_paths = []

    for root, dirs, files in os.walk(module_dir):
        # Skip ignored directories
        dirs[:] = [d for d in dirs if d not in ignore_dirs]
        files = [f for f in files if f not in ignore_files]

        for file in files:
            if file.endswith(".py"):
                script_paths.append(os.path.join(root, file))

    master_dict = {}
    master_class_dict = {}
    master_class_inheritance = {}

    nodes = {}
    node_list = []
    node_idx = 0
    edges = []

    for script_path in script_paths:
        classes, class_dict, methods, orphan_funcs, imports, all_imports, class_inheritance = extract_all(
            script_path, module_dir
        )
        dependency_graph = get_dependency_graph(script_path, module_dir)
        filtered_graph = {}
        for k, v in dependency_graph.items():
            filtered_graph[k] = set()
            for dep in v:
                if dep in classes or dep in methods or dep in orphan_funcs or dep in all_imports:
                    filtered_graph[k].add(dep)

        filtered_graph = {k: v for k, v in filtered_graph.items() if v}

        for class_name in classes:
            nodes[node_idx] = {
                "type": "class",
                "name": class_name,
                "color": NODE_COLORS["class"],
            }
            node_idx += 1
            node_list.append(class_name)

        for method_name in methods:
            nodes[node_idx] = {
                "type": "method",
                "name": method_name,
                "color": NODE_COLORS["method"],
            }
            node_idx += 1
            node_list.append(method_name)

        for func_name in orphan_funcs:
            nodes[node_idx] = {
                "type": "function",
                "name": func_name,
                "color": NODE_COLORS["function"],
            }
            node_idx += 1
            node_list.append(func_name)

        for k, v in class_dict.items():
            for method_name in v:
                if k in node_list and method_name in node_list:
                    edges.append({
                        "from": k,
                        "to": method_name,
                        "type": "defines",
                        "color": EDGE_COLORS["defines"],
                    })
        for k, v in class_inheritance.items():
            for parent_class in v:
                if k in node_list and parent_class in node_list:
                    edges.append({
                        "from": k,
                        "to": parent_class,
                        "type": "inherits",
                        "color": EDGE_COLORS["inherits"],
                    })
        for k, v in filtered_graph.items():
            for dep in v:
                if k in node_list and dep in node_list:
                    edges.append({
                        "from": k,
                        "to": dep,
                        "type": "calls",
                        "color": EDGE_COLORS["calls"],
                    })

    # pprint(node_list)
    # exit()
    G = create_graph(edges, nodes)

    net = Network(
        height="90vh",
        # width="100%",
        # width="15000px",
        directed=True,
        notebook=True,
        layout=False,
        filter_menu=True,
        cdn_resources="in_line",
        # font_color='#10000000',
    )
    # remove orphan nodes (no incoming and no outgoing edges) so they won't appear in the visualization
    orphan_nodes = [n for n in G.nodes() if G.in_degree(n) == 0 and G.out_degree(n) == 0]
    if orphan_nodes:
        G.remove_nodes_from(orphan_nodes)
    G.nodes(data=True)
    net.from_nx(G)
    # net.toggle_physics(True)
    # net.show_buttons(filter_=["physics"])
    options = {
      "physics": {
        "forceAtlas2Based": {
          "gravitationalConstant": -150,
          "springLength": 290
        },
        "minVelocity": 0.75,
        "solver": "forceAtlas2Based"
      }
    }
    net.set_options(json.dumps(options))
    net.show("af_pipeline_mayajaal.html")