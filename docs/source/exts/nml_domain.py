"""
nml_domain.py — Sphinx domain for documenting Fortran namelists.
"""

from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional, Tuple

from docutils.parsers.rst import roles, nodes
from docutils.parsers.rst import directives
from sphinx.domains import Domain
from sphinx.roles import XRefRole
from sphinx.util import logging
from sphinx.util.docutils import SphinxDirective, SphinxRole
from sphinx.util.nodes import make_refnode

logger = logging.getLogger(__name__)


# ============================================================
# Inline roles
# ============================================================

class NmlValueRole(SphinxRole):
    """Render a single value as a code literal."""

    def run(self):
        text = self.text.strip()
        return [nodes.literal(text, text)], []


class NmlValueListRole(SphinxRole):
    """Render a space-separated list of values as literals separated by |."""

    def run(self):
        values = self.text.strip().split()
        out: List[nodes.Node] = []

        for i, val in enumerate(values):
            if i:
                out.append(nodes.Text(" | "))
            out.append(nodes.literal(val, val))

        return out, []


# ============================================================
# Cross-reference roles
# ============================================================

class NmlOptionRole(XRefRole):
    """
    Cross-reference role for options.

    Captures group context at parse time so resolution is reliable.
    """

    def __init__(self) -> None:
        super().__init__(warn_dangling=True)

    def process_link(self, env, refnode, has_explicit_title, title, target):
        refnode["nml:group"] = env.ref_context.get("nml:group")
        return title, target

    def result_nodes(self, document, env, node, is_ref):
        if node.children:
            text = node.astext()
            node.children = [nodes.literal(text, text)]
        return [node], []


class NmlGroupRole(XRefRole):
    """
    Cross-reference role for groups.

    Displays &groupname (literal), resolves against plain name.
    """

    def __init__(self) -> None:
        super().__init__(warn_dangling=True)

    def result_nodes(self, document, env, node, is_ref):
        if node.children:
            group_name = node.astext()
            display = f"&{group_name}"
            node.children = [nodes.literal(display, display)]
        return [node], []


# ============================================================
# Directives
# ============================================================

class NmlGroupDirective(SphinxDirective):
    """
    .. nml:group:: <group_name>
       :no-target:

    Sets the default group (namespace) for subsequent content.

    By default, also registers a hyperlink target for :nml:group:.
    Use :no-target: to suppress target creation.
    """

    required_arguments = 1
    has_content = False

    option_spec = {
        "no-target": directives.flag,
    }

    def run(self) -> List[nodes.Node]:
        env = self.env
        domain: NmlDomain = env.get_domain("nml")  # type: ignore

        group_name = self.arguments[0].strip()

        # Set default namespace for remainder of file
        env.temp_data["nml:group"] = group_name
        env.ref_context["nml:group"] = group_name

        if "no-target" in self.options:
            return []

        target_id = f"nml-group-{group_name}"
        target_node = nodes.target("", "", ids=[target_id])

        source, line = self.get_source_info()
        location = f"{source}:{line}"

        domain.note_group(group_name, env.docname, target_id, location=location)

        return [target_node]


class NmlOptionDirective(SphinxDirective):
    """
    .. nml:option:: <option_name>
       :type:
       :default:

    Defines an option in the current group.
    """

    required_arguments = 1
    has_content = True

    option_spec = {
        "type": directives.unchanged,
        "default": directives.unchanged,
    }

    def run(self) -> List[nodes.Node]:
        env = self.env
        domain: NmlDomain = env.get_domain("nml")  # type: ignore

        option_name = self.arguments[0].strip()

        group = env.temp_data.get("nml:group")
        if not group:
            source, line = self.get_source_info()
            raise ValueError(
                f"Option '{option_name}' defined outside any nml:group "
                f"at {source}:{line}"
            )

        fqname = f"{group}.{option_name}"
        target_id = f"nml-option-{group}-{option_name}"
        target_node = nodes.target("", "", ids=[target_id])

        source, line = self.get_source_info()
        location = f"{source}:{line}"

        domain.note_option(
            fqname,
            env.docname,
            target_id,
            metadata=dict(self.options),
            location=location,
        )

        # Build definition list entry
        dlist = nodes.definition_list()
        item = nodes.definition_list_item()

        term = nodes.term()
        term += nodes.strong("", "", nodes.literal(option_name, option_name))

        if "type" in self.options:
            term += nodes.Text("  (type: ")
            term += nodes.strong(self.options["type"], self.options["type"])
            term += nodes.Text(")")

        if "default" in self.options:
            term += nodes.Text("  default: ")
            dval = self.options["default"]
            term += nodes.strong("", "", nodes.literal(dval, dval))

        definition = nodes.definition()
        self.state.nested_parse(self.content, self.content_offset, definition)

        item += term
        item += definition
        dlist += item

        return [target_node, dlist]


# ============================================================
# Domain
# ============================================================

class NmlDomain(Domain):
    name = "nml"
    label = "Namelist"

    directives = {
        "group": NmlGroupDirective,
        "option": NmlOptionDirective,
    }

    roles = {
        "option": NmlOptionRole(),
        "group": NmlGroupRole(),
    }

    initial_data: Dict[str, Any] = {
        "groups": {},
        "options": {},
    }

    # Incremental rebuild support
    def clear_doc(self, docname: str) -> None:
        groups = self.data.get("groups", {})
        for gname, (dname, _) in list(groups.items()):
            if dname == docname:
                del groups[gname]

        options = self.data.get("options", {})
        for fqname, (dname, _, _) in list(options.items()):
            if dname == docname:
                del options[fqname]

    # Registration
    def note_group(self, name, docname, target_id, location=None):
        groups = self.data["groups"]

        if name in groups:
            existing_docname, _ = groups[name]
            if existing_docname == docname:
                groups[name] = (docname, target_id)
                return

            logger.warning(
                f"Duplicate NML group '{name}' "
                f"(previously defined in {existing_docname})",
                location=location,
            )
            return

        groups[name] = (docname, target_id)

    def note_option(self, fqname, docname, target_id, metadata, location=None):
        options = self.data["options"]

        if fqname in options:
            existing_docname, _, _ = options[fqname]
            if existing_docname == docname:
                options[fqname] = (docname, target_id, metadata)
                return

            logger.warning(
                f"Duplicate NML option '{fqname}' "
                f"(previously defined in {existing_docname})",
                location=location,
            )
            return

        options[fqname] = (docname, target_id, metadata)

    # Xref resolution
    def resolve_xref(
        self, env, fromdocname, builder,
        typ, target, node, contnode
    ) -> Optional[nodes.Node]:

        if typ == "group":
            groups = self.data["groups"]
            if target in groups:
                todocname, target_id = groups[target]
                return make_refnode(
                    builder, fromdocname, todocname,
                    target_id, contnode, target
                )
            return None

        if typ == "option":
            options = self.data["options"]

            # Fully-qualified
            if target in options:
                todocname, target_id, _ = options[target]
                return make_refnode(
                    builder, fromdocname, todocname,
                    target_id, contnode, target
                )

            # Relative to captured group
            group = node.get("nml:group") or env.ref_context.get("nml:group")
            if group:
                fqname = f"{group}.{target}"
                if fqname in options:
                    todocname, target_id, _ = options[fqname]
                    return make_refnode(
                        builder, fromdocname, todocname,
                        target_id, contnode, fqname
                    )

        return None

    def get_objects(self) -> Iterator[Tuple[str, str, str, str, str, int]]:
        for name, (docname, target_id) in self.data["groups"].items():
            yield (name, name, "group", docname, target_id, 1)

        for fqname, (docname, target_id, _) in self.data["options"].items():
            yield (fqname, fqname, "option", docname, target_id, 1)


# ============================================================
# Setup
# ============================================================

def setup(app):
    app.add_domain(NmlDomain)
    app.add_role("nml:value", NmlValueRole())
    app.add_role("nml:valuelist", NmlValueListRole())

    roles.register_generic_role('nml:literal', nodes.literal)

    return {
        "version": "4.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
