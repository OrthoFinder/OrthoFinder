# -*- coding: utf-8 -*-
# #START_LICENSE###########################################################
#
#
# This file is part of the Environment for Tree Exploration program
# (ETE).  http://ete.cgenomics.org
#  
# ETE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ETE is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with ETE.  If not, see <http://www.gnu.org/licenses/>.
#
# 
#                     ABOUT THE ETE PACKAGE
#                     =====================
# 
# ETE is distributed under the GPL copyleft license (2008-2011).  
#
# If you make use of ETE in published work, please cite:
#
# Jaime Huerta-Cepas, Joaquin Dopazo and Toni Gabaldon.
# ETE: a python Environment for Tree Exploration. Jaime BMC
# Bioinformatics 2010,:24doi:10.1186/1471-2105-11-24
#
# Note that extra references to the specific methods implemented in 
# the toolkit are available in the documentation. 
# 
# More info at http://ete.cgenomics.org
#
# 
# #END_LICENSE#############################################################
__VERSION__="ete2-2.2rev1026" 
#START_LICENSE###########################################################
#
# Copyright (C) 2009 by Jaime Huerta Cepas. All rights reserved.
# email: jhcepas@gmail.com
#
# This file is part of the Environment for Tree Exploration program (ETE).
# http://ete.cgenomics.org
#
# ETE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ETE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ETE.  If not, see <http://www.gnu.org/licenses/>.
#
# #END_LICENSE#############################################################
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle
import random
import copy
from collections import deque 
import itertools
import functools
try:
    from rich import print
except ImportError:
    ...

from .newick import read_newick, write_newick

TREEVIEW = False
__all__ = ["Tree", "TreeNode"]

DEFAULT_COMPACT = False
DEFAULT_SHOWINTERNAL = False
DEFAULT_DIST = 1.0
DEFAULT_SUPPORT = 1.0
DEFAULT_NAME = "NoName"

class TreeError(Exception):
    """
    A problem occurred during a TreeNode operation
    """
    def __init__(self, value=''):
        self.value = value
    def __str__(self):
        return repr(self.value)

class TreeNode(object):
    """ 
    TreeNode (Tree) class is used to store a tree structure. A tree
    consists of a collection of TreeNode instances connected in a
    hierarchical way. Trees can be loaded from the New Hampshire Newick
    format (newick).

    :argument newick: Path to the file containing the tree or, alternatively,
       the text string containing the same information.

    :argument 0 format: subnewick format 

      .. table::                                               

          ======  ============================================== 
          FORMAT  DESCRIPTION                                    
          ======  ============================================== 
          0        flexible with support values                  
          1        flexible with internal node names             
          2        all branches + leaf names + internal supports 
          3        all branches + all names                      
          4        leaf branches + leaf names                    
          5        internal and leaf branches + leaf names       
          6        internal branches + leaf names                
          7        leaf branches + all names                     
          8        all names                                     
          9        leaf names                                    
          100      topology only                                 
          ======  ============================================== 

    :returns: a tree node object which represents the base of the tree.

    ** Examples: **

    :: 

        t1 = Tree() # creates an empty tree
        t2 = Tree('(A:1,(B:1,(C:1,D:1):0.5):0.5);')
        t3 = Tree('/home/user/myNewickFile.txt')
    """

    def _get_dist(self):
        return self._dist
    def _set_dist(self, value):
        try:
            self._dist = float(value)
        except ValueError:
            raise

    def _get_support(self):
        return self._support
    def _set_support(self, value):
        try:
            self._support = float(value)
        except ValueError:
            raise

    def _get_up(self):
        return self._up
    def _set_up(self, value):
        if type(value) == type(self) or value is None:
            self._up = value
        else:
            raise ValueError("bad node_up type")

    def _get_children(self):
        return self._children
    def _set_children(self, value):
        if type(value) == list and \
           len(set([type(n)==type(self) for n in value]))<2:
            self._children = value
        else:
            raise ValueError("bad children type")

    def _get_style(self):
        if self._img_style is None:
            self._set_style(None)
           
        return self._img_style

    def _set_style(self, value):
        self.set_style(value)
            
    #: Branch length distance to parent node. Default = 0.0
    img_style = property(fget=_get_style, fset=_set_style)
             
    #: Branch length distance to parent node. Default = 0.0
    dist = property(fget=_get_dist, fset=_set_dist)
    #: Branch support for current node
    support = property(fget=_get_support, fset=_set_support)
    #: Pointer to parent node
    up = property(fget=_get_up, fset=_set_up)
    #: A list of children nodes
    children = property(fget=_get_children, fset=_set_children)

    def _set_face_areas(self, value):
        if isinstance(value, _FaceAreas):
            self._faces = value
        else:
            raise ValueError("[%s] is not a valid FaceAreas instance" %type(value))
        
    def _get_face_areas(self):
        if not hasattr(self, "_faces"):
            self._faces = _FaceAreas()
        return self._faces

    faces = property(fget=_get_face_areas, \
                         fset=_set_face_areas)

    def __init__(self, newick=None, format=0, dist=None, support=None,
                 name=None):
        self._children = []
        self._up = None
        self._dist = DEFAULT_DIST
        self._support = DEFAULT_SUPPORT
        self._img_style = None
        self.features = set([])
        # Add basic features
        self.features.update(["dist", "support", "name"])
        if dist is not None:
            self.dist = dist
        if support is not None:
            self.support = support

        self.name = name if name is not None else DEFAULT_NAME
        
        # Initialize tree
        if newick is not None:
            read_newick(newick, root_node = self, format=format)
           

    def __bool__(self):
        return True
    
    __nonzero__ = __bool__

    def __repr__(self):
        return "Tree node '%s' (%s)" %(self.name, hex(self.__hash__()))

    def __and__(self, value):
        """ This allows to execute tree&'A' to obtain the descendant node
        whose name is A"""
        value=str(value)
        try:
            first_match = next(self.iter_search_nodes(name=value))
            return first_match
        except StopIteration:
            raise ValueError("Node not found")

    def __add__(self, value):
        """ This allows to sum two trees."""
        # Should a make the sum with two copies of the original trees?
        if type(value) == self.__class__:
            new_root = self.__class__()
            new_root.add_child(self)
            new_root.add_child(value)
            return new_root
        else:
            raise ValueError("Invalid node type")

    def __str__(self):
        """ Print tree in newick format. """
        return self.get_ascii(compact=DEFAULT_COMPACT, \
                                show_internal=DEFAULT_SHOWINTERNAL)

    def __contains__(self, item):
        """ Check if item belongs to this node. The 'item' argument must
        be a node instance or its associated name."""
        if isinstance(item, self.__class__):
            return item in set(self.get_descendants())
        elif type(item)==str:
            return item in set([n.name for n in self.traverse()])

    def __len__(self):
        """Node len returns number of children."""
        return len(self.get_leaves())

    def __iter__(self):
        """ Iterator over leaf nodes"""
        return self.iter_leaves()

    def add_feature(self, pr_name, pr_value):
        """ 
        Add or update a node's feature. 
        """
        setattr(self, pr_name, pr_value)
        self.features.add(pr_name)

    def add_features(self, **features):
        """ 
        Add or update several features. """
        for fname, fvalue in features.items():
            setattr(self, fname, fvalue)
            self.features.add(fname)

    def del_feature(self, pr_name):
        """ 
        Permanently deletes a node's feature. 
        """
        if hasattr(self, pr_name):
            delattr(self, pr_name)
            self.features.remove(pr_name)

    # Topology management
    def add_child(self, child=None, name=None, dist=None, support=None):
        """
        Adds a new child to this node. If child node is not suplied
        as an argument, a new node instance will be created.

        :argument None child: the node instance to be added as a child.
        :argument None name: the name that will be given to the child.
        :argument None dist: the distance from the node to the child.
        :argument None support': the support value of child partition.

        :returns: The child node instance

        """
        if child is None:
            child = self.__class__()
        
        if name is not None:
            child.name = name
        if dist is not None:
            child.dist = dist
        if support is not None:
            child.support = support
            
        self.children.append(child)
        child.up = self
        return child

    def remove_child(self, child):
        """ 
        Removes a child from this node (parent and child
        nodes still exit but are no longer connected). 
        """
        try:
            self.children.remove(child)
        except ValueError as e:
            raise TreeError(e)
        else:
            child.up = None
            return child

    def add_sister(self, sister=None, name=None, dist=None):
        """
        Adds a sister to this node. If sister node is not supplied
        as an argument, a new TreeNode instance will be created and
        returned.
        """
        if self.up == None:
            raise TreeError("A parent node is required to add a sister")
        else:
            return self.up.add_child(child=sister, name=name, dist=dist)

    def remove_sister(self, sister=None):
        """
        Removes a node's sister node. It has the same effect as
        **`TreeNode.up.remove_child(sister)`**

        If a sister node is not supplied, the first sister will be deleted
        and returned.

        :argument sister: A node instance

        :return: The node removed
        """
        sisters = self.get_sisters()
        if len(sisters)>0:
            if sister==None:
                sister = sisters.pop(0)
            return self.up.remove_child(sister)

    def delete(self, prevent_nondicotomic=True, preserve_branch_length=False):
        """
        Deletes node from the tree structure. Notice that this method
        makes 'disappear' the node from the tree structure. This means
        that children from the deleted node are transferred to the
        next available parent.

        :param True prevent_nondicotomic: When True (default), delete
        function will be execute recursively to prevent single-child
        nodes.

        :param False preserve_branch_length: If True, branch lengths
        of the deleted nodes are transferred (summed up) to its
        parent's branch, thus keeping original distances among nodes.
                
        **Example:**

        ::

                / C
          root-|
               |        / B
                \--- H |
                        \ A

          > root.delete(H) will produce this structure:

                / C
               |
          root-|--B
               |
                \ A

        """
        parent = self.up
        if parent:
            for ch in self.children:
                if preserve_branch_length:
                    ch.dist += self.dist
                parent.add_child(ch)

            parent.remove_child(self)

        # Avoids parents with only one child
        if prevent_nondicotomic and parent and\
              len(parent.children)<2:
            parent.delete(prevent_nondicotomic=False,
                          preserve_branch_length=preserve_branch_length)


    def detach(self):
        """
        Detachs this node (and all its descendants) from its parent
        and returns the referent to itself.

        Detached node conserves all its structure of descendants, and can
        be attached to another node through the 'add_child' function. This
        mechanism can be seen as a cut and paste.
        """

        if self.up:
            self.up.children.remove(self)
            self.up = None
        return self


    def prune(self, nodes, preserve_branch_length=False):
        """Prunes the topology of a node in order to conserve only a
        selected list of leaf or internal nodes. The minimum number of
        internal nodes (the deepest as possible) are kept to conserve
        the topological relationship among the provided list of nodes.

        :var nodes: a list of node names or node objects that must be kept
        
        :param False preserve_branch_length: If True, branch lengths
        of the deleted nodes are transferred (summed up) to its
        parent's branch, thus keeping original distances among nodes.
        
        **Examples:**

        ::

          t = Tree("(((A:0.1, B:0.01):0.001, C:0.0001):1.0[&&NHX:name=I], (D:0.00001):0.000001[&&NHX:name=J]):2.0[&&NHX:name=root];")
          node_C = t.search_nodes(name="C")[0]
          t.prune(["A","D", node_C])
          print t

        """
        def cmp_nodes(x, y):
            # if several nodes are in the same path of two kept nodes,
            # only one should be maintained. This prioritize internal
            # nodes that are already in the to_keep list and then
            # deeper nodes (closer to the leaves). 
            if x in to_keep:
                if y not in to_keep:
                    return -1
                elif y in to_keep:
                    return 0
            elif n2depth[x] > n2depth[y]:
                return -1
            elif n2depth[x] < n2depth[y]:
                return 1
            else:
                return 0

        to_keep = set(_translate_nodes(self, *nodes))
        start, node2path = self.get_common_ancestor(to_keep, get_path=True)

        # Calculate which kept nodes are visiting the same nodes in
        # their path to the common ancestor.
        n2count = {}
        n2depth = {}
        for seed, path in node2path.items():
            for visited_node in path: 
                if visited_node not in n2depth:
                    depth = visited_node.get_distance(start, topology_only=True)
                    n2depth[visited_node] = depth
                if visited_node is not seed:
                    n2count.setdefault(visited_node, set()).add(seed)

        # if several internal nodes are in the path of exactly the
        # same kept nodes, only one should be maintain. 
        visitors2nodes = {}
        for node, visitors in n2count.items():
            # keep nodes connection at least two other nodes
            if len(visitors)>1: 
                visitor_key = frozenset(visitors)
                visitors2nodes.setdefault(visitor_key, set()).add(node)
        for visitors, nodes in visitors2nodes.items():
            s = sorted(nodes, key=functools.cmp_to_key(cmp_nodes))
            to_keep.add(s[0])

        # Detach unvisited branches
        if start is not self:
            start.detach()
            for n in self.get_children(): 
                n.detach()
            for n in start.get_children():
                self.add_child(child=n)

        for n in [self]+self.get_descendants():
            if n not in to_keep: 
                n.delete(prevent_nondicotomic=False,
                         preserve_branch_length=preserve_branch_length)

    def swap_children(self):
        """
        Swaps current children order.
        """
        if len(self.children)>1:
            self.children.reverse()

    # def prune_OLD(self, nodes):
    #     """
    #     Prunes the topology of this node in order to conserve only a
    #     selected list of leaf or internal nodes. The algorithm deletes
    #     nodes until getting a consistent topology with a subset of
    #     nodes. Topology relationships among kept nodes is maintained.
    #  
    #     :var nodes: a list of node names or node objects that must be kept
    #  
    #     **Examples:**
    #  
    #     ::
    #  
    #       t = Tree("(((A:0.1, B:0.01):0.001, C:0.0001):1.0[&&NHX:name=I], (D:0.00001):0.000001[&&NHX:name=J]):2.0[&&NHX:name=root];")
    #       node_C = t.search_nodes(name="C")[0]
    #       t.prune(["A","D", node_C])
    #       print t
    #     """
    #    
    #     to_keep = set(_translate_nodes(self, *nodes))
    #     to_detach = []
    #     for node in self.traverse("postorder"):
    #         for c in node.children:
    #             if c in to_keep:
    #                 to_keep.add(node)
    #                 break
    #         if node not in to_keep:
    #             to_detach.append(node)
    #             for c in node.children:
    #                 to_detach.remove(c)
    #     for node in to_detach:
    #         node.detach()
    #     for node in to_keep:
    #         if len(node.children) == 1:
    #             node.delete()
    #     if len(self.children)==1 and self.children[0] not in to_keep:
    #         self.children[0].delete()

    # #####################
    # Tree traversing 
    # #####################


    def get_children(self):
        """ 
        Returns an independent list of node's children. 
        """
        return [ch for ch in self.children]

    def get_sisters(self):
        """ 
        Returns an indepent list of sister nodes. 
        """
        if self.up!=None:
            return [ch for ch in self.up.children if ch!=self]
        else:
            return []

    def iter_leaves(self, is_leaf_fn=None):
        """ 
        Returns an iterator over the leaves under this node. 

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.traverse(strategy="preorder", is_leaf_fn=is_leaf_fn):
            if not is_leaf_fn:
                if n.is_leaf():
                    yield n
            else:
                if is_leaf_fn(n):
                    yield n
                    
    def get_leaves(self, is_leaf_fn=None):
        """
        Returns the list of terminal nodes (leaves) under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [n for n in self.iter_leaves(is_leaf_fn=is_leaf_fn)]

    def iter_leaf_names(self, is_leaf_fn=None):
        """ 
        Returns an iterator over the leaf names under this node. 

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.iter_leaves(is_leaf_fn=is_leaf_fn):
            yield n.name

    def get_leaf_names(self, is_leaf_fn=None):
        """
        Returns the list of terminal node names under the current
        node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [name for name in self.iter_leaf_names(is_leaf_fn=is_leaf_fn)]

    def iter_descendants(self, strategy="levelorder", is_leaf_fn=None):
        """ 
        Returns an iterator over all descendant nodes. 

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.traverse(strategy=strategy, is_leaf_fn=is_leaf_fn):
            if n is not self:
                yield n

    def get_descendants(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns a list of all (leaves and internal) descendant nodes.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [n for n in self.iter_descendants(strategy=strategy, \
                                                 is_leaf_fn=is_leaf_fn)]

    def delete_traverse(self):
        for node in self.traverse("postorder"):
            if not node.is_root():
                node.up.remove_child(node)
            # node.delete()
            del node

    def traverse(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns an iterator to traverse the tree structure under this
        node.
         
        :argument "levelorder" strategy: set the way in which tree
           will be traversed. Possible values are: "preorder" (first
           parent and then children) 'postorder' (first children and
           the parent) and "levelorder" (nodes are visited in order
           from root to leaves)

        :argument None is_leaf_fn: If supplied, ``is_leaf_fn``
           function will be used to interrogate nodes about if they
           are terminal or internal. ``is_leaf_fn`` function should
           receive a node instance as first argument and return True
           or False. Use this argument to traverse a tree by
           dynamically collapsing internal nodes matching
           ``is_leaf_fn``.
        """
        if strategy=="preorder":
            return self._iter_descendants_preorder(is_leaf_fn=is_leaf_fn)
        elif strategy=="levelorder":
            return self._iter_descendants_levelorder(is_leaf_fn=is_leaf_fn)
        elif strategy=="postorder":
            return self._iter_descendants_postorder(is_leaf_fn=is_leaf_fn)
            
    def _iter_descendants_postorder_recursive(self, is_leaf_fn=None):
        """
        Iterate over all desdecendant nodes.
        """
        if not is_leaf_fn or not is_leaf_fn(self):
            for ch in self.children:
                for node in ch._iter_descendants_postorder(is_leaf_fn=is_leaf_fn):
                    yield node
        yield self
    
    def iter_prepostorder(self, is_leaf_fn=None):
        """
        Iterate over all nodes in a tree yielding every node in both
        pre and post order. Each iteration returns a postorder flag
        (True if node is being visited in postorder) and a node
        instance.
        """
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                yield (False, node)
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
            else:
                #POSTORDER ACTIONS
                yield (True, node)

    def _iter_descendants_postorder(self, is_leaf_fn=None):
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
                else:
                    yield node
            else:
                #POSTORDER ACTIONS
                yield node
        
    def _iter_descendants_levelorder(self, is_leaf_fn=None):
        """ 
        Iterate over all desdecendant nodes. 
        """
        tovisit = deque([self])
        while len(tovisit)>0:
            node = tovisit.popleft()
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                tovisit.extend(node.children)

    def _iter_descendants_preorder(self, is_leaf_fn=None):
        """
        Iterator over all descendant nodes. 
        """
        to_visit = deque()
        node = self
        while node is not None:
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                to_visit.extendleft(reversed(node.children))
            try:
                node = to_visit.popleft()
            except:
                node = None

    def iter_ancestors(self):
        '''versionadded: 2.2
        
        Iterates over the list of all ancestor nodes from current node
        to the current tree root.

        '''
        node = self
        while node.up is not None:
            yield node.up
            node = node.up

    def get_ancestors(self):
        '''versionadded: 2.2

        Returns the list of all ancestor nodes from current node to
        the current tree root.

        '''
        return [n for n in self.iter_ancestors()]
            
    def describe(self):
        """ 
        Prints general information about this node and its
        connections.
        """
        if len(self.get_tree_root().children)==2:
            rooting = "Yes"
        elif len(self.get_tree_root().children)>2:
            rooting = "No"
        else:
            rooting = "Unknown"
        max_node, max_dist = self.get_farthest_leaf()
        cached_content = self.get_cached_content()
        print("Number of leaf nodes:\t%d" % len(cached_content[self]))
        print("Number of internal nodes:\t%d" % len(cached_content))
        print("Rooted:\t%s" %rooting)
        print("Most distant node:\t%s" %max_node.name)
        print("Max. distance:\t%f" %max_dist)
        
    def write(self, features=None, outfile=None, format=0, is_leaf_fn=None,
              format_root_node=False):
        """ 
        Returns the newick representation of current node. Several
        arguments control the way in which extra data is shown for
        every node:

        :argument features: a list of feature names to be exported
          using the Extended Newick Format (i.e. features=["name",
          "dist"]). Use an empty list to export all available features
          in each node (features=[])

        :argument outfile: writes the output to a given file

        :argument format: defines the newick standard used to encode the
          tree. See tutorial for details.

        :argument False format_root_node: If True, it allows features
          and branch information from root node to be exported as a
          part of the newick text string. For newick compatibility
          reasons, this is False by default.

        :argument is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
          
        **Example:**

        ::

             t.get_newick(features=["species","name"], format=1)

        """

        nw = write_newick(self, features = features, format=format,
                          is_leaf_fn=is_leaf_fn,
                          format_root_node=format_root_node)
        if outfile is not None:
            open(outfile, "w").write(nw)
        else:
            return nw

    def get_tree_root(self):
        """ 
        Returns the absolute root node of current tree structure. 
        """
        root = self
        while root.up is not None:
            root = root.up
        return root

    def get_common_ancestor(self, *target_nodes, **kargs):
        """ 
        Returns the first common ancestor between this node and a given
        list of 'target_nodes'.

        **Examples:**

        ::

          t = tree.Tree("(((A:0.1, B:0.01):0.001, C:0.0001):1.0[&&NHX:name=common], (D:0.00001):0.000001):2.0[&&NHX:name=root];")
          A = t.get_descendants_by_name("A")[0]
          C = t.get_descendants_by_name("C")[0]
          common =  A.get_common_ancestor(C)
          print common.name

        """
        
        get_path = kargs.get("get_path", False)

        if len(target_nodes) == 1 and type(target_nodes[0]) \
                in set([set, tuple, list, frozenset]):
            target_nodes = target_nodes[0]

        # Convert node names into node instances
        target_nodes = _translate_nodes(self, *target_nodes)

        # If only one node is provided, use self as the second target
        if type(target_nodes) != list:
            target_nodes = [target_nodes, self]
        elif len(target_nodes)==1:
            target_nodes = tree_nodes.append(self)

        n2path = {}
        reference = []
        ref_node = None
        for n in target_nodes:
            current = n
            while current: 
                n2path.setdefault(n, set()).add(current)
                if not ref_node:
                    reference.append(current)
                current = current.up
            if not ref_node:
                ref_node = n

        common = None
        for n in reference:
            broken = False
            for node, path in n2path.items():
                if node is not ref_node and n not in path:
                    broken = True
                    break

            if not broken: 
                common = n
                break
        if not common: 
            raise ValueError("Nodes are not connected!")

        if get_path:
            return common, n2path
        else:
            return common

    def iter_search_nodes(self, **conditions):
        """ 
        Search nodes in an interative way. Matches are being yield as
        they are being found. This avoids to scan the full tree
        topology before returning the first matches. Useful when
        dealing with huge trees.
        """
        
        for n in self.traverse():
            conditions_passed = 0
            for key, value in conditions.items():
                if hasattr(n, key) and getattr(n, key) == value:
                    conditions_passed +=1
            if conditions_passed == len(conditions):
                yield n

    def search_nodes(self, **conditions):
        """
        Returns the list of nodes matching a given set of conditions.

        **Example:**

        ::

          tree.search_nodes(dist=0.0, name="human")

        """
        matching_nodes = []
        for n in self.iter_search_nodes(**conditions):
            matching_nodes.append(n)
        return matching_nodes

    def get_leaves_by_name(self,name):
        """ 
        Returns a list of leaf nodes matching a given name. 
        """
        return self.search_nodes(name=name, children=[])

    def is_leaf(self):
        """ 
        Return True if current node is a leaf.
        """
        return len(self.children) == 0

    def is_root(self):
        """ 
        Returns True if current node has no parent
        """
        if self.up is None:
            return True
        else:
            return False

    # ###########################
    # Distance related functions
    # ###########################
    def get_distance(self, target, target2=None, topology_only=False):
        """
        Returns the distance between two nodes. If only one target is
        specified, it returns the distance bewtween the target and the
        current node.
        
        :argument target: a node within the same tree structure.

        :argument target2: a node within the same tree structure. If
          not specified, current node is used as target2.

        :argument False topology_only: If set to True, distance will
          refer to the number of nodes between target and target2.

        :returns: branch length distance between target and
          target2. If topology_only flag is True, returns the number
          of nodes between target and target2.
 
        """

        if target2 is None:
            target2 = self
            root = self.get_tree_root()
        else:
            # is target node under current node?
            root = self

        target, target2 = _translate_nodes(root, target, target2)
        ancestor = root.get_common_ancestor(target, target2)
        if ancestor is None:
            raise TreeError("Nodes are not connected")

        dist = 0.0
        for n in [target2, target]:
            current = n
            while current != ancestor:
                if topology_only:
                    if  current!=target:
                        dist += 1
                else:
                    dist += current.dist
                current = current.up
        return dist

    def get_farthest_node(self, topology_only=False):
        """
        Returns the node's farthest descendant or ancestor node, and the
        distance to it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the farthest node referred to the
          current node and the distance to it.

        """
        # Init fasthest node to current farthest leaf
        farthest_node,farthest_dist = self.get_farthest_leaf(topology_only=topology_only)
        prev    = self
        if topology_only:
            cdist = 0
        else:
            cdist = prev.dist
        current = prev.up
        while current is not None:
            for ch in current.children:
                if ch != prev:
                    if not ch.is_leaf():
                        fnode, fdist = ch.get_farthest_leaf(topology_only=topology_only)
                    else:
                        fnode = ch
                        fdist = 0
                    if topology_only:
                        fdist += 1.0
                    else:
                        fdist += ch.dist
                    if cdist+fdist > farthest_dist:
                        farthest_dist = cdist + fdist
                        farthest_node = fnode
            prev = current
            if topology_only:
                cdist += 1
            else:
                cdist  += prev.dist
            current = prev.up
        return farthest_node, farthest_dist

    def get_farthest_leaf(self, topology_only=False):
        """
        Returns node's farthest descendant node (which is always a leaf), and the
        distance to it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the farthest leaf referred to the
          current node and the distance to it.
        """
        max_dist = 0.0
        max_node = None
        if self.is_leaf():
            return self, 0.0
        else:
            for ch in self.children:
                node, d = ch.get_farthest_leaf(topology_only=topology_only)
                if topology_only:
                    d += 1.0
                else:
                    d += ch.dist
                if d>=max_dist:
                    max_dist = d
                    max_node = node
            return max_node, max_dist

    def get_closest_leaf(self, topology_only=False):
        """Returns node's closest descendant leaf and the distance to
        it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the closest leaf referred to the
          current node and the distance to it.

        """
        min_dist = None
        min_node = None
        if self.is_leaf():
            return self, 0.0
        else:
            for ch in self.children:
                node, d = ch.get_closest_leaf(topology_only=topology_only)
                if topology_only:
                    d += 1.0
                else:
                    d += ch.dist
                if min_dist is None or d<min_dist:
                    min_dist = d
                    min_node = node
            return min_node, min_dist


            
    def get_midpoint_outgroup(self):
        """
        Returns the node that divides the current tree into two distance-balanced
        partitions.
        """
        # Gets the farthest node to the current root
        root = self.get_tree_root()
        nA , r2A_dist = root.get_farthest_leaf()
        nB , A2B_dist = nA.get_farthest_node()

        outgroup = nA
        middist  = A2B_dist / 2.0
        cdist = 0
        current = nA
        while current is not None:
            cdist += current.dist
            if cdist > (middist): # Deja de subir cuando se pasa del maximo
                break
            else:
                current = current.up
        return current
        
    def populate(self, size, names_library=None, reuse_names=False,
                 random_branches=False, branch_range=(0,1),
                 support_range=(0,1)): 
        """
        Generates a random topology by populating current node.

        :argument None names_library: If provided, names library
          (list, set, dict, etc.) will be used to name nodes.

        :argument False reuse_names: If True, node names will not be
          necessarily unique, which makes the process a bit more
          efficient.

        :argument False random_branches: If True, branch distances and support
          values will be randomized.
        
        :argument (0,1) branch_range: If random_branches is True, this
        range of values will be used to generate random distances.

        :argument (0,1) support_range: If random_branches is True,
        this range of values will be used to generate random branch
        support values.

        """
        NewNode = self.__class__

        if len(self.children) > 1: 
            connector = NewNode()
            for ch in self.get_children():
                ch.detach()
                connector.add_child(child = ch)
            root = NewNode()
            self.add_child(child = connector)
            self.add_child(child = root)
        else:
            root = self

        next = deque([root])
        for i in range(size-1):
            if random.randint(0, 1):
                p = next.pop()
            else:
                p = next.popleft()

            c1 = p.add_child()
            c2 = p.add_child()
            next.extend([c1, c2])
            if random_branches:
                c1.dist = random.uniform(*branch_range)
                c2.dist = random.uniform(*branch_range)
                c1.support = random.uniform(*branch_range)
                c2.support = random.uniform(*branch_range)
            else:
                c1.dist = 1.0
                c2.dist = 1.0
                c1.support = 1.0
                c2.support = 1.0

        # next contains leaf nodes
        charset =  "abcdefghijklmnopqrstuvwxyz"
        if names_library:
            names_library = deque(names_library)
        else:
            avail_names = itertools.combinations_with_replacement(charset, 10)
        for n in next:
            if names_library:
                if reuse_names: 
                    tname = random.sample(names_library, 1)[0]
                else:
                    tname = names_library.pop()
            else:
                tname = ''.join(next(avail_names))
            n.name = tname
            

    def set_outgroup(self, outgroup):
        """
        Sets a descendant node as the outgroup of a tree.  This function
        can be used to root a tree or even an internal node.

        :argument outgroup: a node instance within the same tree
          structure that will be used as a basal node. 

        """

        outgroup = _translate_nodes(self, outgroup)

        if self == outgroup:
            raise ValueError("Cannot set myself as outgroup")

        parent_outgroup = outgroup.up

        # Detects (sub)tree root
        n = outgroup
        while n.up is not self:
            n = n.up

        # If outgroup is a child from root, but with more than one
        # sister nodes, creates a new node to group them

        self.children.remove(n)
        if len(self.children) != 1:
            down_branch_connector = self.__class__()
            down_branch_connector.dist = 0.0
            down_branch_connector.support = n.support
            for ch in self.get_children():
                down_branch_connector.children.append(ch)
                ch.up = down_branch_connector
                self.children.remove(ch)
        else:
            down_branch_connector = self.children[0]

        # Connects down branch to myself or to outgroup
        quien_va_ser_padre = parent_outgroup
        if quien_va_ser_padre is not self:
            # Parent-child swapping
            quien_va_ser_hijo = quien_va_ser_padre.up
            quien_fue_padre = None
            buffered_dist = quien_va_ser_padre.dist
            buffered_support = quien_va_ser_padre.support
            
            while quien_va_ser_hijo is not self:
                quien_va_ser_padre.children.append(quien_va_ser_hijo)
                quien_va_ser_hijo.children.remove(quien_va_ser_padre)

                buffered_dist2 = quien_va_ser_hijo.dist
                buffered_support2 = quien_va_ser_hijo.support
                quien_va_ser_hijo.dist = buffered_dist
                quien_va_ser_hijo.support = buffered_support
                buffered_dist = buffered_dist2
                buffered_support = buffered_support2

                quien_va_ser_padre.up = quien_fue_padre
                quien_fue_padre = quien_va_ser_padre

                quien_va_ser_padre = quien_va_ser_hijo
                quien_va_ser_hijo = quien_va_ser_padre.up

            quien_va_ser_padre.children.append(down_branch_connector)
            down_branch_connector.up = quien_va_ser_padre
            quien_va_ser_padre.up = quien_fue_padre

            down_branch_connector.dist += buffered_dist
            outgroup2 = parent_outgroup
            parent_outgroup.children.remove(outgroup)
            outgroup2.dist = 0
           
        else:
            outgroup2 = down_branch_connector

        outgroup.up = self
        outgroup2.up = self
        # outgroup is always the first children. Some function my
        # trust on this fact, so do no change this.
        self.children = [outgroup,outgroup2]
        middist = (outgroup2.dist + outgroup.dist)/2
        outgroup.dist = middist
        outgroup2.dist = middist
        outgroup2.support = outgroup.support

    def unroot(self):
        """ 
        Unroots current node. This function is expected to be used on
        the absolute tree root node, but it can be also be applied to
        any other internal node. It will convert a split into a
        multifurcation.
        """
        # if is rooted
        if not self.is_root():
            sys.stderr.write("Warning. You are unrooting an internal node!!!\n")
        if len(self.children)==2:
            if not self.children[0].is_leaf():
                self.children[0].delete()
            elif not self.children[1].is_leaf():
                self.children[1].delete()
            else:
                raise TreeError("Cannot unroot a tree with only two leaves")

    def copy(self, method="cpickle"):
        """.. versionadded: 2.1

        Returns a copy of the current node.

        :var cpickle method: Protocol used to copy the node
        structure. The following values are accepted:

           - "newick": Tree topology, node names, branch lengths and
             branch support values will be copied by as represented in
             the newick string (copy by newick string serialisation).
        
           - "newick-extended": Tree topology and all node features
             will be copied based on the extended newick format
             representation. Only node features will be copied, thus
             excluding other node attributes. As this method is also
             based on newick serialisation, features will be converted
             into text strings when making the copy. 

           - "cpickle": The whole node structure and its content is
             cloned based on cPickle object serialisation (slower, but
             recommended for full tree copying)
        
           - "deepcopy": The whole node structure and its content is
             copied based on the standard "copy" Python functionality
             (this is the slowest method but it allows to copy complex
             objects even if attributes point to lambda functions,
             etc.)

        """
        if method=="newick":
            new_node = self.__class__(self.write(features=["name"], format_root_node=True))
        elif method=="newick-extended":
            self.write(features=[], format_root_node=True)
            new_node = self.__class__(self.write(features=[]))
        elif method == "deepcopy":
            parent = self.up
            self.up = None
            new_node = copy.deepcopy(self)
            self.up = parent
        elif method == "cpickle":
            parent = self.up
            self.up = None
            new_node = pickle.loads(pickle.dumps(self, 2))
            self.up = parent
        else:
            raise ValuerError("Invalid copy method")
            
        return new_node
        
    def _asciiArt(self, char1='-', show_internal=True, compact=False, attributes=None):
        """
        Returns the ASCII representation of the tree.

        Code based on the PyCogent GPL project.
        """
        if not attributes:
            attributes = ["name"]
        node_name = ', '.join(map(str, [getattr(self, v) for v in attributes if hasattr(self, v)]))
        
        LEN = max(3, len(node_name) if not self.children or show_internal else 3)
        PAD = ' ' * LEN
        PA = ' ' * (LEN-1)
        if not self.is_leaf():
            mids = []
            result = []
            for c in self.children:
                if len(self.children) == 1:
                    char2 = '/'
                elif c is self.children[0]:
                    char2 = '/'
                elif c is self.children[-1]:
                    char2 = '\\'
                else:
                    char2 = '-'
                (clines, mid) = c._asciiArt(char2, show_internal, compact, attributes)
                mids.append(mid+len(result))
                result.extend(clines)
                if not compact:
                    result.append('')
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = [PAD] * (lo+1) + [PA+'|'] * (hi-lo-1) + [PAD] * (end-hi)
            mid = int((lo + hi) / 2)
            prefixes[mid] = char1 + '-'*(LEN-2) + prefixes[mid][-1]
            result = [p+l for (p,l) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + node_name + stem[len(node_name)+1:]
            return (result, mid)
        else:
            return ([char1 + '-' + node_name], 0)

    def get_ascii(self, show_internal=True, compact=False, attributes=None):
        """
        Returns a string containing an ascii drawing of the tree.

        :argument show_internal: includes internal edge names.
        :argument compact: use exactly one line per tip.

        :param attributes: A list of node attributes to shown in the
            ASCII representation.
        
        """
        (lines, mid) = self._asciiArt(show_internal=show_internal,
                                      compact=compact, attributes=attributes)
        return '\n'+'\n'.join(lines)


    def ladderize(self, direction=0):
        """ 
        .. versionadded: 2.1 

        Sort the branches of a given tree (swapping children nodes)
        according to the size of each partition.

        ::

           t =  Tree("(f,((d, ((a,b),c)),e));")

           print t

           #            
           #      /-f
           #     |
           #     |          /-d
           # ----|         |
           #     |     /---|          /-a
           #     |    |    |     /---|
           #     |    |     \---|     \-b
           #      \---|         |
           #          |          \-c
           #          |
           #           \-e

           t.ladderize()
           print t

           #      /-f
           # ----|
           #     |     /-e
           #      \---|
           #          |     /-d
           #           \---|
           #               |     /-c
           #                \---|
           #                    |     /-a
           #                     \---|
           #                          \-b

        """
       
        if not self.is_leaf():
            n2s = {}
            for n in self.get_children():
                s = n.ladderize(direction=direction)
                n2s[n] = s

            self.children.sort(lambda x,y: cmp(n2s[x], n2s[y]))
            if direction == 1:
                self.children.reverse()
            size = sum(n2s.values())
        else:
            size = 1

        return size

    def sort_descendants(self):
        """ 
        .. versionadded: 2.1 

        This function sort the branches of a given tree by
        considerening node names. After the tree is sorted, nodes are
        labeled using ascendent numbers.  This can be used to ensure
        that nodes in a tree with the same node names are always
        labeled in the same way. Note that if duplicated names are
        present, extra criteria should be added to sort nodes.

        Unique id is stored as a node._nid attribute
       
        """

        node2content = self.get_cached_content()
        def sort_by_content(x, y):
            return cmp(str(sorted([i.name for i in node2content[x]])),
                       str(sorted([i.name for i in node2content[y]])))

        for n in self.traverse():
            if not n.is_leaf():
                n.children.sort(sort_by_content)
        return node2content

    # def get_cached_content(self, store_attr=None,  _store=None):
    #     """ 
    #     .. versionadded: 2.2
       
    #     Returns a dictionary pointing to the preloaded content of each
    #     internal node under this tree. Such a dictionary is intended
    #     to work as a cache for operations that require many traversal
    #     operations.

    #     :param None store_attr: Specifies the node attribute that
    #     should be cached (i.e. name, distance, etc.). When none, the
    #     whole node instance is cached.

    #     :param _store: (internal use)

    #     """
    #     if _store is None:
    #         _store = {}
            
    #     for ch in self.children:
    #         ch.get_cached_content(store_attr=store_attr, _store=_store)

    #     if self.children:
    #         val = set()
    #         for ch in self.children:
    #             val.update(_store[ch])
    #         _store[self] = val
    #     else:
    #         if store_attr is None:
    #             val = self
    #         else:
    #             val = getattr(self, store_attr)
    #         _store[self] = set([val])
    #     return _store

    
    # def get_cached_content_non_recursive(self, store_attr=None,  _store=None):
    def get_cached_content(self, store_attr=None,  _store=None):
        """ 
        Returns a dictionary pointing to the preloaded content of each
        internal node under this tree. Such a dictionary is intended
        to work as a cache for operations that require many traversal
        operations.

        :param None store_attr: Specifies the node attribute that
        should be cached (i.e. name, distance, etc.). When none, the
        whole node instance is cached.

        :param _store: (internal use)

        """
        if _store is None:
            _store = {}

        for n in self.traverse("postorder"):
            if n.children:
                val = set()
                for ch in n.children:
                    val.update(_store[ch])
                _store[n] = val
            else:
                if store_attr is None:
                    val = n
                else:
                    val = getattr(n, store_attr)
                _store[n] = set([val])
        return _store

       
    def robinson_foulds(self, t2, attr_t1="name", attr_t2="name"):
        """
        .. versionadded: 2.2
        
        Returns the Robinson-Foulds symmetric distance between current
        tree and a different tree instance.
     
        :param t2: target tree
        
        :param name attr_t1: Compare trees using a custom node
                              attribute as a node name.
        
        :param name attr_t2: Compare trees using a custom node
                              attribute as a node name in target tree.

        :returns: (symmetric distance, total partitions, common node
         names, partitions in current tree, partitions in target tree)
           
        """
        
        t1 = self
        t1content = t1.get_cached_content()
        t2content = t2.get_cached_content()
        target_names = set([getattr(_n, attr_t1) for _n in t1content[t1]])
        ref_names = set([getattr(_n, attr_t2) for _n in t2content[t2]])
        common_names = target_names & ref_names
        if len(common_names) < 2:
            raise ValueError("Trees share less than 2 nodes")

        r1 = set([",".join(sorted([getattr(_c, attr_t1) for _c in cont
                                   if getattr(_c, attr_t1) in common_names]))
                  for cont in t1content.values() if len(cont)>1])
        r2 = set([",".join(sorted([getattr(_c, attr_t2) for _c in cont
                                   if getattr(_c, attr_t2) in common_names]))
                  for cont in t2content.values() if len(cont)>1])
        r1.discard("")
        r2.discard("")              
        inters = r1.intersection(r2)
        if len(r1) == len(r2):
                rf = (len(r1) - len(inters)) * 2
        else :
                rf = (len(r1) - len(inters)) + (len(r2) - len(inters))
        max_parts = len(r1) + len(r2)
        return rf, max_parts, common_names, r1, r2

    def get_partitions(self):
        """ 
        .. versionadded: 2.1
        
        It returns the set of all possible partitions under a
        node. Note that current implementation is quite inefficient
        when used in very large trees.

        t = Tree("((a, b), e);")
        partitions = t.get_partitions()

        # Will return: 
        # a,b,e
        # a,e
        # b,e
        # a,b
        # e
        # b
        # a
        """
        all_leaves = frozenset(self.get_leaf_names())
        all_partitions = set([all_leaves])
        for n in self.iter_descendants():
            p1 = frozenset(n.get_leaf_names())
            p2 = frozenset(all_leaves - p1)
            all_partitions.add(p1)
            all_partitions.add(p2)
        return all_partitions

    def convert_to_ultrametric(self, tree_length, strategy="balanced"):
        """
        .. versionadded: 2.1 

        Converts a tree to ultrametric topology (all leaves must have
        the same distance to root). Note that, for visual inspection
        of ultrametric trees, node.img_style["size"] should be set to
        0.
        """

        # pre-calculate how many splits remain under each node
        node2max_depth = {}
        for node in self.traverse("postorder"):
            if not node.is_leaf():
                max_depth = max([node2max_depth[c] for c in node.children]) + 1
                node2max_depth[node] = max_depth
            else:
                node2max_depth[node] = 1
        node2dist = {self: 0.0}
        tree_length = float(tree_length)
        step = tree_length / node2max_depth[self]
        for node in self.iter_descendants("levelorder"):
            if strategy == "balanced":
                node.dist = (tree_length - node2dist[node.up]) / node2max_depth[node]
                node2dist[node] =  node.dist + node2dist[node.up]
            elif strategy == "fixed":
                if not node.is_leaf():
                    node.dist = step
                else:
                    node.dist = tree_length - ((node2dist[node.up]) * step)
                node2dist[node] = node2dist[node.up] + 1
            node.dist = node.dist

    def check_monophyly(self, values, target_attr, ignore_missing=False):
        """
        Returns True if a given target attribute is monophyletic under
        this node for the provided set of values. 

        If not all values are represented in the current tree
        structure, a ValueError exception will be raised to warn that
        strict monophyly could never be reached (this behaviour can be
        avoided by enabling the `ignore_missing` flag.

        :param values: a set of values for which monophyly is
            expected.
            
        :param target_attr: node attribute being used to check
            monophyly (i.e. species for species trees, names for gene
            family trees, or any custom feature present in the tree).

        :param False ignore_missing: Avoid raising an Exception when
            missing attributes are found. 
        """
        if type(values) != set:
            values = set(values)

        # This is the only time I traverse the tree, then I use cached
        # leaf content
        n2leaves = self.get_cached_content()
        # Locate leaves matching requested attribute values
        targets = [leaf for leaf in n2leaves[self]
                   if getattr(leaf, target_attr) in values]
       
        # Raise an error if requested attribute values are not even present
        if not ignore_missing:
            missing_values = values - set([getattr(n, target_attr) for n
                                           in targets])
            if missing_values: 
                raise ValueError("Expected '%s' value(s) not found: %s" %(
                    target_attr, ','.join(missing_values)))
                
        # Check monophyly with get_common_ancestor. Note that this
        # step does not require traversing the tree again because
        # targets are node instances instead of node names, and
        # get_common_ancestor function is smart enough to detect it
        # and avoid unnecessary traversing.
        common = self.get_common_ancestor(targets)
        observed = n2leaves[common]
        foreign_leaves = [leaf for leaf in observed
                          if getattr(leaf, target_attr) not in values]
        if not foreign_leaves:
            return True, "monophyletic"
        else:
            # if the requested attribute is not monophyletic in this
            # node, let's differentiate between poly and paraphyly. 
            poly_common = self.get_common_ancestor(foreign_leaves)
            # if the common ancestor of all foreign leaves is self
            # contained, we have a paraphyly. Otherwise, polyphyly. 
            polyphyletic = [leaf for leaf in poly_common if
                            getattr(leaf, target_attr) in values]
            if polyphyletic:
                return False, "polyphyletic"
            else:
                return False, "paraphyletic"

            
    def get_monophyletic(self, values, target_attr):
        """
        .. versionadded:: 2.2

        Returns a list of nodes matching the provided monophyly
        criteria. For a node to be considered a match, all
        `target_attr` values within and node, and exclusively them,
        should be grouped.
        
        :param values: a set of values for which monophyly is
            expected.
            
        :param target_attr: node attribute being used to check
            monophyly (i.e. species for species trees, names for gene
            family trees).
           
        """

        if type(values) != set:
            values = set(values)

        n2values = self.get_cached_content(store_attr=target_attr)
      
        is_monophyletic = lambda node: n2values[node] == values
        for match in self.iter_leaves(is_leaf_fn=is_monophyletic):
            if is_monophyletic(match):
                yield match

            
    def resolve_polytomy(self, default_dist=0.0, default_support=0.0,
                         recursive=True):
        """
        .. versionadded: 2.2
        
        Resolve all polytomies under current node by creating an
        arbitrary dicotomic structure among the affected nodes. This
        function randomly modifies current tree topology and should
        only be used for compatibility reasons (i.e. programs
        rejecting multifurcated node in the newick representation).

        :param 0.0 default_dist: artificial branch distance of new
            nodes.
                                  
        :param 0.0 default_support: artificial branch support of new
            nodes.
                                   
        :param True recursive: Resolve any polytomy under this
             node. When False, only current node will be checked and fixed.
        """
        
        
        def _resolve(node):
            if len(node.children) > 2:
                children = list(node.children)
                node.children = []
                next_node = root = node
                for i in range(len(children)-2):
                    next_node = next_node.add_child()
                    next_node.dist = default_dist
                    next_node.support = default_support

                next_node = root
                for ch in children:
                    next_node.add_child(ch)
                    if ch != children[-2]:
                        next_node = next_node.children[0]
        target = [self]
        if recursive: 
            target.extend([n for n in self.get_descendants()])
        for n in target:
            _resolve(n)


def _translate_nodes(root, *nodes):
    name2node = dict([ [n, None] for n in nodes if type(n) is str])
    for n in root.traverse():
        if n.name in name2node:
            if name2node[n.name] is not None:
                raise ValueError("Ambiguous node name: "+str(n.name))
            else:
                name2node[n.name] = n

    if None in name2node.values():
        notfound = [key for key, value in name2node.items() if value is None]
        raise ValueError("Node names not found: "+str(notfound))

    valid_nodes = []
    for n in nodes: 
        if type(n) is not str:
            if type(n) is not root.__class__ :
                raise ValueError("Invalid target node: "+str(n))
            else:
                valid_nodes.append(n)
            
    valid_nodes.extend(name2node.values())
    if len(valid_nodes) == 1:
        return valid_nodes[0]
    else:
        return valid_nodes


# Alias
#: .. currentmodule:: ete2
Tree = TreeNode




    
