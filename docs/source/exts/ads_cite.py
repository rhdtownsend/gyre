import pickle
import ads
import re
from docutils.parsers.rst import nodes

## Role functions for citing ADS articles

def build_cite(rawtext, ref, lineno, inliner, options, template):

    if not ref in ads_data:
        msg = inliner.reporter.error(
            'Reference {:s} not found'.format(ref), line=lineno)
        prb = inliner.problematic(rawtext, rawtext, msg)
        return [prb], [msg]

    year_str = re.findall(r'[0-9]+[A-Za-z]*', ref)

    if len(year_str) > 0:
        year = year_str[0]
    else:
        year = ads_data[ref].year
        
    if len(ads_data[ref].author) == 1:
        author = format(ads_data[ref].author[0].split(',')[0])
    elif len(ads_data[ref].author) == 2:
        author = '{:s} & {:s}'.format(ads_data[ref].author[0].split(',')[0],
                                      ads_data[ref].author[1].split(',')[0])
    else:
        author = '{:s} et al.'.format(ads_data[ref].author[0].split(',')[0])

    citation = template.format(author, year)
    url = 'https://ui.adsabs.harvard.edu/abs/{:s}/abstract'.format(ads_data[ref].bibcode)

    node = nodes.reference(rawtext, citation, refuri=url, **options)

    return [node], []


def ads_citet(role, rawtext, text, lineno, inliner,
              options={}, content=[]):

    return build_cite(rawtext, text, lineno, inliner, options, '{0:s} ({1:s})')


def ads_citep(role, rawtext, text, lineno, inliner,
              options={}, content=[]):

    return build_cite(rawtext, text, lineno, inliner, options, '({0:s}, {1:s})')


def ads_citealt(role, rawtext, text, lineno, inliner,
                options={}, content=[]):

    return build_cite(rawtext, text, lineno, inliner, options, '{0:s} {1:s}')


def ads_citealp(role, rawtext, text, lineno, inliner,
                options={}, content=[]):

    return build_cite(rawtext, text, lineno, inliner, options, '{0:s}, {1:s}')


def ads_citeauthor(role, rawtext, text, lineno, inliner,
                   options={}, content=[]):

    return build_cite(rawtext, text, lineno, inliner, options, '{0:s}')


def ads_citeyear(role, rawtext, text, lineno, inliner,
                options={}, content=[]):

    return build_cite(rawtext, text, lineno, inliner, options, '{1:s}')


def setup(app):

    # Read the data

    global ads_data

    with open('{:s}/ads_refs.dat'.format(app.srcdir), 'rb') as f:
        ads_data = pickle.load(f).copy()

    # Set up roles

    app.add_role('ads_citet', ads_citet)
    app.add_role('ads_citep', ads_citep)
    app.add_role('ads_citealt', ads_citealt)
    app.add_role('ads_citealp', ads_citealp)
    app.add_role('ads_citeauthor', ads_citeauthor)
    app.add_role('ads_citeyear', ads_citeyear)
    
    return {
        'version': '0.1',
    }

