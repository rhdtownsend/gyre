"""
ofile_domain.py — Sphinx domain for documenting output files.
"""

from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional, Tuple

#from docutils import nodes
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

# class OfileValueRole(SphinxRole):
#     """Render a single value as a code literal."""

#     def run(self):
#         text = self.text.strip()
#         return [nodes.literal(text, text)], []


# class OfileValueListRole(SphinxRole):
#     """Render a space-separated list of values as literals separated by |."""

#     def run(self):
#         values = self.text.strip().split()
#         out: List[nodes.Node] = []

#         for i, val in enumerate(values):
#             if i:
#                 out.append(nodes.Text(" | "))
#             out.append(nodes.literal(val, val))

#         return out, []


# ============================================================
# Cross-reference roles
# ============================================================

class OfileFieldRole(XRefRole):
    """
    Cross-reference role for fields.

    Captures filetype context at parse time so resolution is reliable.
    """

    def __init__(self) -> None:
        super().__init__(warn_dangling=True)

    def process_link(self, env, refnode, has_explicit_title, title, target):
        refnode["ofile:filetype"] = env.ref_context.get("ofile:filetype")
        return title, target

    def result_nodes(self, document, env, node, is_ref):
        if node.children:
            text = node.astext()
            node.children = [nodes.literal(text, text)]
        return [node], []


class OfileFiletypeRole(XRefRole):
    """
    Cross-reference role for filetype.
    """

    def __init__(self) -> None:
        super().__init__(warn_dangling=True)

    def result_nodes(self, document, env, node, is_ref):
        if node.children:
            filetype_name = node.astext()
            display = filetype_name
            node.children = [nodes.literal(display, display)]
        return [node], []


# ============================================================
# Directives
# ============================================================

class OfileFiletypeDirective(SphinxDirective):
    """
    .. ofile:filetype:: <filetype_name>
       :no-target:

    Sets the default filetype (namespace) for subsequent content.

    By default, also registers a hyperlink target for :ofile:filetype:.
    Use :no-target: to suppress target creation.
    """

    required_arguments = 1
    has_content = False

    option_spec = {
        "no-target": directives.flag,
    }

    def run(self) -> List[nodes.Node]:
        env = self.env
        domain: OfileDomain = env.get_domain("ofile")  # type: ignore

        filetype_name = self.arguments[0].strip()

        # Set default namespace for remainder of file
        env.temp_data["ofile:filetype"] = filetype_name
        env.ref_context["ofile:filetype"] = filetype_name

        if "no-target" in self.options:
            return []

        target_id = f"ofile-filetype-{filetype_name}"
        target_node = nodes.target("", "", ids=[target_id])

        source, line = self.get_source_info()
        location = f"{source}:{line}"

        domain.note_filetype(filetype_name, env.docname, target_id, location=location)

        return [target_node]


class OfileFieldDirective(SphinxDirective):
    """
    .. ofile:field:: <field_name>
       :type:
       :dim:

    Defines a field in the current filetype.
    """

    required_arguments = 1
    has_content = True

    option_spec = {
        "type": directives.unchanged,
        "dim": directives.unchanged,
        "units": directives.unchanged,
    }

    def run(self) -> List[nodes.Node]:

        msgs: List[nodes.system_message] = []

        env = self.env
        domain: OfileDomain = env.get_domain("ofile")  # type: ignore

        field_name = self.arguments[0].strip()

        if "type" not in self.options:
            raise self.error(
                "The :type: option is required for .. nml:field::"
            )

        filetype = env.temp_data.get("ofile:filetype")
        if not filetype:
            source, line = self.get_source_info()
            raise ValueError(
                f"Field '{field_name}' defined outside any ofile:filetype "
                f"at {source}:{line}"
            )

        fqname = f"{filetype}.{field_name}"
        target_id = f"ofile-field-{filetype}-{field_name}"
        target_node = nodes.target("", "", ids=[target_id])

        source, line = self.get_source_info()
        location = f"{source}:{line}"

        domain.note_field(
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
        term += nodes.strong("", "", nodes.literal(field_name, field_name))

        term += [nodes.Text(" ("),
                 nodes.emphasis("", "", nodes.Text("type:")),
                 nodes.Text(f' {self.options["type"]}')]

        if "dim" in self.options:
            term += [nodes.Text(", "),
                     nodes.emphasis("", "", nodes.Text("dimension:")),
                     nodes.Text(" ")]
            dval = self.options["dim"]
            dim_nodes, dim_msgs = self.state.inline_text(dval, self.lineno)

            term += dim_nodes

            msgs.extend(dim_msgs)

        if "units" in self.options:
            term += [nodes.Text(", "),
                     nodes.emphasis("", "", nodes.Text("units:")),
                     nodes.Text(" ")]
            dval = self.options["units"]
            units_nodes, units_msgs = self.state.inline_text(dval, self.lineno)

            term += units_nodes

            msgs.extend(units_msgs)

        term += nodes.Text(")")

        definition = nodes.definition()
        self.state.nested_parse(self.content, self.content_offset, definition)

        item += term
        item += definition
        dlist += item

        return [target_node, dlist]


# ============================================================
# Domain
# ============================================================

class OfileDomain(Domain):
    name = "ofile"
    label = "Outfile"

    directives = {
        "filetype": OfileFiletypeDirective,
        "field": OfileFieldDirective,
    }

    roles = {
        "field": OfileFieldRole(),
        "filetype": OfileFiletypeRole(),
    }

    initial_data: Dict[str, Any] = {
        "filetypes": {},
        "fields": {},
    }

    # Incremental rebuild support
    def clear_doc(self, docname: str) -> None:
        filetype = self.data.get("filetype", {})
        for name, (dname, _) in list(filetype.items()):
            if dname == docname:
                del filetype[name]

        fields = self.data.get("fields", {})
        for name, (dname, _, _) in list(fields.items()):
            if dname == docname:
                del fields[name]

    # Registration
    def note_filetype(self, name, docname, target_id, location=None):
        filetypes = self.data["filetypes"]

        if name in filetypes:
            existing_docname, _ = filetypes[name]
            if existing_docname == docname:
                filetypes[name] = (docname, target_id)
                return

            logger.warning(
                f"Duplicate OFILE filetype '{name}' "
                f"(previously defined in {existing_docname})",
                location=location,
            )
            return

        filetypes[name] = (docname, target_id)

    def note_field(self, fqname, docname, target_id, metadata, location=None):
        fields = self.data["fields"]

        if fqname in fields:
            existing_docname, _, _ = fields[fqname]
            if existing_docname == docname:
                fields[fqname] = (docname, target_id, metadata)
                return

            logger.warning(
                f"Duplicate OFILE field '{fqname}' "
                f"(previously defined in {existing_docname})",
                location=location,
            )
            return

        fields[fqname] = (docname, target_id, metadata)

    # Xref resolution
    def resolve_xref(
        self, env, fromdocname, builder,
        typ, target, node, contnode
    ) -> Optional[nodes.Node]:

        if typ == "filetype":
            filetypes = self.data["filetypes"]
            if target in filetypes:
                todocname, target_id = filetypes[target]
                return make_refnode(
                    builder, fromdocname, todocname,
                    target_id, contnode, target
                )
            return None

        if typ == "field":
            fields = self.data["fields"]

            # Fully-qualified
            if target in fields:
                todocname, target_id, _ = fields[target]
                return make_refnode(
                    builder, fromdocname, todocname,
                    target_id, contnode, target
                )

            # Relative to captured group
            filetype = node.get("ofile:filetype") or env.ref_context.get("ofile:filetype")
            if filetype:
                fqname = f"{filetype}.{target}"
                if fqname in fields:
                    todocname, target_id, _ = fields[fqname]
                    return make_refnode(
                        builder, fromdocname, todocname,
                        target_id, contnode, fqname
                    )

        return None

    def get_objects(self) -> Iterator[Tuple[str, str, str, str, str, int]]:
        for name, (docname, target_id) in self.data["filetypes"].items():
            yield (name, name, "filetype", docname, target_id, 1)

        for fqname, (docname, target_id, _) in self.data["fields"].items():
            yield (fqname, fqname, "field", docname, target_id, 1)


# ============================================================
# Setup
# ============================================================

def setup(app):
    app.add_domain(OfileDomain)
    #app.add_role("ofile:value", OfileValueRole())
    #app.add_role("ofile:valuelist", OfileValueListRole())

    roles.register_generic_role('ofile:literal', nodes.literal)

    return {
        "version": "4.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
