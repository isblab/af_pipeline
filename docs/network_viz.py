import ast
import os
from pprint import pprint
from collections import defaultdict
import networkx as nx
from pyvis.network import Network
import json
import argparse

NODE_COLORS = {
    "class": "#FFDD00",
    "method": "#00A3FF",
    "function": "#00FF19",
    # "import": "#FF00F6",
    "unknown": "#D3D3D3",
}

EDGE_COLORS = {
    "defines": "#E3E3E3",
    "calls": "#274c77",
    "inherits": "#ff595e",
    "has_a": "#9d4edd",
}

REPO_LINK = "https://github.com/isblab/af_pipeline/tree/main/af_pipeline"

COMMON_IGNORES = {
    "ignore_dirs": [
        "docs",
        "tests",
        "utils",
    ],
    "ignore_files": [],
}

MODULE_SPECIFIC_IGNORES = {
    "af_pipeline": {
        "ignore_dirs": [],
        "ignore_files": [],
    },
    "af_input": {
        "ignore_dirs": [
            # "af_input",
            "constants",
            "pae_to_domains",
            "parser",
            "rank_predictions",
            "rigid_bodies",
            "tools",
        ],
        "ignore_files": [
            "_initialize.py",
        ],
    },
    "pae_to_domains": {
        "ignore_dirs": [
            "af_input",
            "constants",
            # "pae_to_domains",
            "parser",
            "rank_predictions",
            "rigid_bodies",
            "tools",
        ],
        "ignore_files": [
            "_initialize.py",
        ],
    },
    "parser": {
        "ignore_dirs": [
            "af_input",
            "constants",
            "pae_to_domains",
            # "parser",
            "rank_predictions",
            "rigid_bodies",
            "tools",
        ],
        "ignore_files": [
            "_initialize.py",
        ],
    },
    "rank_predictions": {
        "ignore_dirs": [
            "af_input",
            "constants",
            "pae_to_domains",
            "parser",
            # "rank_predictions",
            "rigid_bodies",
            "tools",
        ],
        "ignore_files": [
            "_initialize.py",
        ],
    },
    "rigid_bodies": {
        "ignore_dirs": [
            "af_input",
            "constants",
            "pae_to_domains",
            "parser",
            "rank_predictions",
            # "rigid_bodies",
            "tools",
        ],
        "ignore_files": [
            "_initialize.py",
        ],
    },
    "interaction": {
        "ignore_dirs": [
            "af_input",
            "constants",
            "pae_to_domains",
            "parser",
            "rank_predictions",
            "rigid_bodies",
            # "tools",
        ],
    },
    "tools": {
        "ignore_dirs": [
            "af_input",
            "constants",
            "pae_to_domains",
            "parser",
            "rank_predictions",
            "rigid_bodies",
            # "tools",
        ],
        "ignore_files": [
            "_initialize.py",
        ],
    },
}

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
    class_compositions = defaultdict(list)


    curr_module = os.path.splitext(script_path.replace(project_dir, ""))[0].replace(os.path.sep, ".")
    class_linenos = {}
    method_linenos = {}
    orphan_funcs_linenos = {}

    if module_name is None:
        module_name = project_dir.split(os.path.sep)[-1]

    curr_module = f"{module_name}{curr_module}"

    with open(script_path, "r") as file:
        tree = ast.parse(file.read())

    for node in ast.walk(tree):

        if isinstance(node, ast.ClassDef): # class definition
            all_classes.append(f"{curr_module}.{node.name}")
            class_linenos[f"{curr_module}.{node.name}"] = node.lineno
            class_dict[f"{curr_module}.{node.name}"] = []
            class_inheritance[f"{curr_module}.{node.name}"] = []
            # print(ast.dump(node, indent=4))
            # exit()
            for item in node.body:
                if isinstance(item, ast.FunctionDef): # method definition
                    all_methods.append(f"{curr_module}.{node.name}.{item.name}")
                    method_linenos[f"{curr_module}.{node.name}.{item.name}"] = item.lineno
                    class_dict[f"{curr_module}.{node.name}"].append(
                        f"{curr_module}.{node.name}.{item.name}"
                    )
            for base in node.bases: # inheritance
                if isinstance(base, ast.Name):
                    class_inheritance[f"{curr_module}.{node.name}"].append(base.id)

        elif isinstance(node, ast.FunctionDef): # function definition (includes methods)
            if node.name.startswith("__") is False:
                orphan_funcs.append(f"{curr_module}.{node.name}")
                orphan_funcs_linenos[f"{curr_module}.{node.name}"] = node.lineno

            elif node.name.startswith("__init__"):

                # find the enclosing class for this __init__ method
                parent_class = None
                for parent in ast.walk(tree):
                    if isinstance(parent, ast.ClassDef):
                        if node in parent.body:
                            parent_class = parent
                            break

                if parent_class is None:
                    continue
                full_class_name = f"{curr_module}.{parent_class.name}"

                # look for instance attribute assignments to capture composition
                for item in node.body:
                    if isinstance(item, ast.Assign):
                        for target in item.targets:
                            if isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name) and target.value.id == 'self':
                                # e.g., self.attr = ClassName()
                                if isinstance(item.value, ast.Call):
                                    if isinstance(item.value.func, ast.Name):
                                        class_name = item.value.func.id
                                        look_up = True
                                        # get full name from imports or defined classes
                                        full_name = f"{curr_module}.{class_name}"
                                        all_classes_now = list(class_dict.keys())
                                        all_classes_now = [
                                            c for c in all_classes_now if c.startswith(curr_module + ".")
                                        ]
                                        if full_name in all_classes_now:
                                            class_name = full_name
                                            look_up = False
                                        if look_up:
                                            # look up for class_name among import_dict for curr_module
                                            all_imports_now = import_dict[curr_module]
                                            for imp in all_imports_now:
                                                if imp.endswith(f".{class_name}"):
                                                    class_name = imp
                                                    break
                                        class_compositions[f"{full_class_name}"].append(class_name)
                                    elif isinstance(item.value.func, ast.Attribute):
                                        parts = []
                                        f = item.value.func
                                        while isinstance(f, ast.Attribute):
                                            parts.append(f.attr)
                                            f = f.value
                                        if isinstance(f, ast.Name):
                                            parts.append(f.id)
                                        parts.reverse()
                                        class_name = ".".join(parts)
                                        look_up = True
                                        # get full name from imports or defined classes
                                        full_name = f"{curr_module}.{class_name}"
                                        all_classes_now = list(class_dict.keys())
                                        all_classes_now = [
                                            c for c in all_classes_now if c.startswith(curr_module + ".")
                                        ]
                                        if full_name in all_classes_now:
                                            class_name = full_name
                                            look_up = False
                                        if look_up:
                                            # look up for class_name among import_dict for curr_module
                                            all_imports_now = import_dict[curr_module]
                                            for imp in all_imports_now:
                                                if imp.endswith(f".{class_name}"):
                                                    class_name = imp
                                                    break
                                        class_compositions[f"{full_class_name}"].append(class_name)

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
    orphan_funcs_linenos = {k: v for k, v in orphan_funcs_linenos.items() if k in orphan_funcs}
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

    return (
        all_classes,
        class_dict,
        all_methods,
        orphan_funcs,
        import_dict,
        all_imports,
        class_inheritance,
        class_compositions,
        class_linenos,
        method_linenos,
        orphan_funcs_linenos,
    )


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

        # relation tracking (e.g., composition / "has_a")
        # self.relations[owner_class]['has_a'] = set(of component class full names)
        self.relations = defaultdict(lambda: defaultdict(set))

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

                        # If this assignment occurs inside __init__, record a "has_a" (composition) relation.
                        # Also add to dependencies so it appears in downstream filtered graphs if desired.
                        if self.current_scope and self.current_scope[-1].endswith('.__init__'):
                            if target_class_full:
                                self.relations[parent_class]['has_a'].add(target_class_full)
                                if parent_class not in self.dependencies:
                                    self.dependencies[parent_class] = set()
                                # add the class as a dependency target (semantic relation captured separately in self.relations)
                                if target_class_full != parent_class:
                                    self.dependencies[parent_class].add(target_class_full)
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

def get_submodule_script_paths(submodule, module_dir):

    ignore_dirs = MODULE_SPECIFIC_IGNORES.get(submodule, {}).get("ignore_dirs", [])
    ignore_files = MODULE_SPECIFIC_IGNORES.get(submodule, {}).get("ignore_files", [])

    ignore_dirs += COMMON_IGNORES.get("ignore_dirs", [])
    ignore_files += COMMON_IGNORES.get("ignore_files", [])

    ignore_dirs = [os.path.join(module_dir, i_dir) for i_dir in ignore_dirs]

    script_paths = []

    for root, dirs, files in os.walk(module_dir):

        files = [f for f in files if f not in ignore_files]
        full_file_paths = [os.path.join(root, f) for f in files]

        if root in ignore_dirs:
            continue

        for file_path in full_file_paths:
            if file_path.endswith(".py"):
                script_paths.append(file_path)

    return script_paths

def get_nodes_and_edges(script_paths_dict, module_dir):

    node_edge_dict = {
        submodule: {
            "nodes": {},
            "edges": [],
        } for submodule in script_paths_dict.keys()
    }

    for submodule, script_paths in script_paths_dict.items():

        nodes = {}
        edges = []
        node_idx = 0
        node_list = []
        all_methods = []
        all_classes = []
        filtered_graph = {}

        for script_path in script_paths + script_paths:

            (
                classes,
                class_dict,
                methods,
                orphan_funcs,
                imports,
                all_imports,
                class_inheritance,
                class_compositions,
                class_linenos,
                method_linenos,
                orphan_funcs_linenos,
            ) = extract_all(
                script_path, module_dir
            )

            dependency_graph = get_dependency_graph(script_path, module_dir)
            class_compositions = dict(class_compositions)

            for k, v in dependency_graph.items():
                filtered_graph[k] = set() if k not in filtered_graph else filtered_graph[k]
                for dep in v:
                    if "save_rb_assessment" in dep and "save_rigid_bodies" in k:
                        print("Considering call edge:", k, "->", dep)
                    if (
                        dep in classes + methods + orphan_funcs + all_imports + all_methods
                    ) or (
                        dep in class_compositions
                    ):
                        filtered_graph[k].add(dep)

            filtered_graph = {k: v for k, v in filtered_graph.items() if v}

            for class_name in classes:
                nodes[node_idx] = {
                    "type": "class",
                    "name": class_name,
                    "color": NODE_COLORS["class"],
                    "lineno": class_linenos.get(class_name, None),
                }
                node_idx += 1
                node_list.append(class_name)

            for method_name in methods:
                all_methods.append(method_name) if method_name not in all_methods else None
                if method_name in node_list:
                    continue
                nodes[node_idx] = {
                    "type": "method",
                    "name": method_name,
                    "color": NODE_COLORS["method"],
                    "lineno": method_linenos.get(method_name, None),
                }
                node_idx += 1
                node_list.append(method_name)

            for func_name in orphan_funcs:
                if func_name in node_list:
                    continue
                nodes[node_idx] = {
                    "type": "function",
                    "name": func_name,
                    "color": NODE_COLORS["function"],
                    "lineno": orphan_funcs_linenos.get(func_name, None),
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
                            "font": {"color": EDGE_COLORS["defines"]},
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

            for k, v in class_compositions.items():
                for dep in v:
                    if k in classes+all_imports and dep in classes+all_imports:
                        last_dep = dep.split(".")[-1]
                        # check capitalization of the first letter to decide if it's a class
                        is_class = last_dep[0].isupper()
                        if is_class:
                            edges.append({
                                "from": k,
                                "to": dep,
                                "type": "has_a",
                                "color": EDGE_COLORS["has_a"],
                            })

            for k, v in filtered_graph.items():
                if "RigidBodies.save_rigid_bodies" in k:
                    print("Processing call edges for:", k)
                    pprint(f"    {v}", )
                for dep in v:
                    if "save_rb_assessment" in dep:# and "save_rigid_bodies" in k:
                        print("###########call edge:", k, "->", dep)
                    if k in node_list + all_methods and dep in node_list + all_methods:
                        edges.append({
                            "from": k,
                            "to": dep,
                            "type": "calls",
                            "color": EDGE_COLORS["calls"],
                        })

        node_edge_dict[submodule]["nodes"] = nodes
        node_edge_dict[submodule]["edges"] = edges

    return node_edge_dict


def generate_network_visualization(edges, nodes, submodule_name, output_dir, remove_orphans=True):

    G = create_graph(edges, nodes)

    net = Network(
        height="90vh",
        width="98vw",
        directed=True,
        notebook=False,
        layout=False,
        # select_menu=True,
        # filter_menu=True,
        cdn_resources="in_line",
        # font_color='#10000000',
    )

    # remove orphan nodes
    if remove_orphans:
        orphan_nodes = [
            n for n in G.nodes() if G.in_degree(n) == 0 and G.out_degree(n) == 0
        ]
        if orphan_nodes:
            G.remove_nodes_from(orphan_nodes)

    # edit node label to show only the basename
    for node_id, data in G.nodes(data=True):
        link_to_add = ""
        node_data = nodes[node_id]

        if data["type"] == "class" or data["type"] == "function":
            parent_script_path = "/".join(data["label"].split(".")[:-1]) + ".py"
            link_to_add = f"{REPO_LINK}/{parent_script_path}#L{node_data.get('lineno', '')}"

        elif data["type"] == "method":
            parent_full_class = ".".join(data["label"].split(".")[:-1])
            parent_script_path = "/".join(parent_full_class.split(".")[:-1]) + ".py"
            link_to_add = f"{REPO_LINK}/{parent_script_path}#L{node_data.get('lineno', '')}"

        data["title"] = (
            f"<a href='{link_to_add}' target='_blank'>{data["label"]}</a>"
        )
        data["label"] = data["label"].split(".")[-1]


    # G.nodes(data=True)
    net.from_nx(G)

    options = {
      "physics": {
        "forceAtlas2Based": {
          "theta": 0.6,
          "gravitationalConstant": -99,
          "springLength": 170,
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
        "tooltipDelay": 300,
      },
      "layout": {
        "hierarchical": {
        "enabled": False,
        "direction": "UD",
        "sortMethod": "directed",
        "shakeTowards": "leaves",
        "nodeSpacing": 150,
        "levelSeparation": 250,
        "treeSpacing": 60,
        }
      },
    }

    net.set_options(json.dumps(options))
    net.set_template_dir(
        "./docs/template",
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
    with open(f"{output_dir}/network_{submodule_name}.html", "w") as f:
        f.write(net_html)
    # net.write_html(f"{output_dir}/network_{submodule_name}.html")

if __name__ == "__main__":

    args = argparse.ArgumentParser()

    args.add_argument(
        "-d",
        "--module_dir",
        type=str,
        # default="/home/omg/Projects/af_pipeline/af_pipeline",
        required=True,
        help="Path to the module directory.",
    )

    args.add_argument(
        "-m",
        "--module_name",
        type=str,
        default="af_pipeline",
        help="Name of the module.",
    )

    args.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default="./docs/network_viz",
        help="Directory to save the network visualizations.",
    )

    args = args.parse_args()

    module_dir = os.path.abspath(args.module_dir)
    module_name = args.module_name

    networks_to_generate = [
        "af_pipeline",
        "af_input",
        # "pae_to_domains",
        "interaction",
        "parser",
        "rank_predictions",
        "rigid_bodies",
        "tools",
    ]

    script_paths_dict = {}

    for submodule in networks_to_generate:
        script_paths_dict[submodule] = get_submodule_script_paths(submodule, module_dir)

    # pprint(script_paths_dict)

    node_edge_dict = get_nodes_and_edges(script_paths_dict, module_dir)

    for submodule, data in node_edge_dict.items():
        edges = data["edges"]
        nodes = data["nodes"]
        generate_network_visualization(
            edges,
            nodes,
            submodule,
            args.output_dir,
            remove_orphans=True,
        )
