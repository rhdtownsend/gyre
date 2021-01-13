from docutils.parsers.rst import roles, nodes

def nml_group_role (name, rawtext, text, lineno, inliner, options={}, content=[]):

    """Format a namelist group name.

    Returns 2 part tuple containing list of nodes to insert into the
    document and a list of system messages.  Both are allowed to be
    empty.

    :param name: The role name used in the document.
    :param rawtext: The entire markup snippet, with role.
    :param text: The text marked with the role.
    :param lineno: The line number where rawtext appears in the input.
    :param inliner: The inliner instance that called us.
    :param options: Directive options for customization.
    :param content: The directive content for customization.
    """

    node = nodes.literal(rawtext, '&'+text)
    return [node], []

def setup (app):

    """Install the plugin.

    :param app: Sphinx application context.
    """

    app.add_role('nml_g', nml_group_role)

    roles.register_generic_role('nml_n', nodes.literal)
    roles.register_generic_role('nml_v', nodes.literal)
    roles.register_generic_role('nml_nv', nodes.literal)

    return
