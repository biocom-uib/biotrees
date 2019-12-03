import requests
import feedparser
from os import path
from urllib.parse import urljoin, quote as urlquote
from xml.etree import ElementTree
from networkx import DiGraph

from biotrees.phylotree import digraph_to_phylotree


"""
Stores host and base path information about the TreeBase API
"""
class TreeBaseClient(object):
    def __init__(self, host, base_path):
        self.host = host
        self.base_path = base_path

    @property
    def base_url(self):
        return urljoin(self.host, self.base_path)


"""
Create a new `TreeBaseClient`
"""
def new_client(host='https://purl.org', base_path='/phylo/treebase/phylows/'):
    return TreeBaseClient(host, base_path)


"""
Retrieve and parse an RDF file describing the search results for trees satisfying a query string.
:param client: `TreeBaseClient` instance
:param query: PhyloWS search query string <https://github.com/TreeBASE/treebase/wiki/API#searching>. Defaults to all Species trees.
:param cache_file: File to read/save to skip re-downloading the same search.
:return: Parsed feed (with `feedparser`)
"""
def search_trees(client, query="tb.kind.tree=Species", cache_file=None):
    url = urljoin(client.base_url, 'tree/find')

    if cache_file is not None and path.isfile(cache_file):
        return feedparser.parse(cache_file)

    else:
        text = requests.get(url, params={'query': query, 'format': 'rss1', 'recordSchema': 'tree'}).text

        with open(cache_file, 'w+') as f:
            f.write(text)

        return feedparser.parse(text)

"""
Specialization of `search_trees` that looks for all trees of kind `Species`.
:param client: `TreeBaseClient` instance
:param cache_file: File to read/save to skip re-downloading the same search.
:return: Parsed feed (with `feedparser`)
"""
def search_all_species_trees(client, cache_file=None):
    return search_trees(client, "tb.kind.tree=Species", cache_file=cache_file)


"""
Specialization of `search_trees` that looks for all trees of kind `Gene`.
:param client: `TreeBaseClient` instance
:param cache_file: File to read/save to skip re-downloading the same search.
:return: Parsed feed (with `feedparser`)
"""
def search_all_gene_trees(client, cache_file=None):
    return search_trees(client, "tb.kind.tree=Gene", cache_file=cache_file)


"""
Get the TreeBase tree IDs from the RDF feed obtained via `search_trees`.
:param feed: The result of `search_trees`
:return: Generator of tree IDs (strings)
"""
def get_feed_tree_ids(feed):
    for item in feed['items']:
        yield item['summary']


"""
Fetch the raw data of a tree from TreeBase
:param client: `TreeBaseClient` instance
:param tree_id: TreeBase tree ID (string)
:param format: Format in which the data should be retrieved
:return: A `(ok, resp)` tuple containing whether the request was successful and the response object if it failed or the response text if succeeded.
"""
def fetch_tree(client, tree_id, format='nexml'):
    url = urljoin(client.base_url, 'tree/')
    url = urljoin(url, urlquote(tree_id))

    response = requests.get(url, params={'format': format})

    if 'accessviolation' in response.url:
        return False, response
    else:
        return True, response.text


"""
Apply `fetch_tree` to each TreeBase tree ID. Yield only successful requests and call `on_error` on error responses.
:param client: `TreeBaseClient` instance
:param tree_id: Iterable of TreeBase tree IDs (strings)
:param format: Format in which the data should be retrieved
:param on_error: Function taking a `tree_id` and a `requests.Response` to be called on failed requests.
:return: Generator of `(tree_id, data)` tuples where data is the raw tree data in the specified format.
"""
def fetch_trees(client, tree_ids, format='nexml', on_error=None):
    for tree_id in tree_ids:
        ok, resp = fetch_tree(client, tree_id, format=format)

        if ok:
            yield tree_id, resp
        elif on_error is not None:
            on_error(tree_id, resp)


"""
Build a `networkx.DiGraph` object from a `<tree>` NeXML element.
"""
def _nexml_tree_to_digraph(tree):
    g = DiGraph()

    for node in tree:
        if node.tag != 'node' and not node.tag.endswith('}node'):
            continue
        g.add_node(node.get('id'))

    for edge in tree:
        if edge.tag != 'edge' and not edge.tag.endswith('}edge'):
            continue
        g.add_edge(edge.get('source'), edge.get('target'))

    return g


"""
Parse an NeXML string and yield a `networkx.DiGraph` for each `<tree>` in it.
:param nexml: NeXML string
:return: Generator of `networkx.DiGraph`s
"""
def nexml_digraphs(nexml):
    for trees in ElementTree.fromstring(nexml):
        if trees.tag != 'trees' and not trees.tag.endswith('}trees'):
            continue

        for tree in trees:
            if tree.tag != 'tree' and not tree.tag.endswith('}tree'):
                continue

            yield tree.get('id'), _nexml_tree_to_digraph(tree)


"""
Convenience function combining `fetch_trees`, `nexml_digraphs` and `digraph_to_phylotree`.
:param client: `TreeBaseClient` instance
:param tree_ids: Iterable of TreeBase tree IDs
:param on_error: Function taking a `tree_id` and a `requests.Response` to be called on individual failed tree requests.
:return: Generator of `(tree_id, phylotree)` pairs.
"""
def fetch_phylotrees(client, tree_ids, on_error=None):
    for tree_id, nexml in fetch_trees(client, tree_ids, format='nexml', on_error=on_error):
        for tree_id, g in nexml_digraphs(nexml):
            yield tree_id, digraph_to_phylotree(g)


if __name__ == '__main__':
    # example printing the shapes of all species trees in TreeBase (takes a long time)

    client = new_client()
    trees_feed = search_trees(client, cache_file='treebase_index.rdf')

    failed = []
    succeeded = []

    from biotrees.phylotree import phylotree_to_shape

    def log_tree(t=None):
        global trees_feed
        global failed
        global succeeded

        if t is not None:
            print(phylotree_to_shape(t))
        else:
            print('{} succeeded {} failed out of {} trees'.format(len(succeeded), len(failed), len(trees_feed['items'])), end='')


    def tree_failed(tree_id, response):
        global failed
        failed.append(tree_id)
        log_tree()


    for tree_id, tree in fetch_phylotrees(client, get_feed_tree_ids(trees_feed), on_error=tree_failed):
        succeeded.append(tree_id)
        log_tree(tree)
