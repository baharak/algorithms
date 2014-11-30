/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package longestPalindrome;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import longestPalindrome.SuffixTree.Node;

public class LongestPalindrome {
	
	/**
	 * Finds the longest substring of s that is a palindrome in O(s.length()).
	 * If there are several of equal length, returns the leftmost one.
	 * Note: due to restriction of SuffixTree library, s can 
	 * only contains lower case letters and digits (no space).
	 * 
	 * @param s a string consisting of lowercase letters and digits only
	 * @return a substring of s (could be empty if s is empty)
	 */
	public static String findIn(String s) {
		int n = s.length();
		int maxPaliLen = 0;
		int[] fromAndTo = new int[]{0, 0};
		LongestCommonPrefix lcp = new LongestCommonPrefix(
				s,
				new StringBuilder(s).reverse().toString());

		// O(n) preprocessing
		List<int[]> lcpQueries = new ArrayList<int[]>(2 * n);
		for (int i = 0; i < n; ++i) {
			lcpQueries.add(new int[]{i, n - i - 1}); // precompute for odd
			lcpQueries.add(new int[]{i, n - i}); // precompute for even
		}
		LongestCommonPrefix.Result lcpResult = lcp.preprocess(lcpQueries); 

		for (int i = 0; i < n; ++i) {
			// check for odd-length palindrome centered at pos i in O(1)
			int k = lcpResult.get(lcpQueries.get(2 * i));
			int paliLen = 2 * k - 1;
			if (paliLen > maxPaliLen) {
				maxPaliLen = paliLen;
				fromAndTo[0] = i - k + 1;
				fromAndTo[1] = i + k;
			}
			// check for even-length palindrome centered at pos i in O(1)
			k = lcpResult.get(lcpQueries.get(2 * i + 1));
			paliLen = 2 * k;
			if (paliLen > maxPaliLen) {
				maxPaliLen = paliLen;
				fromAndTo[0] = i - k;
				fromAndTo[1] = i + k;
			}
		}
		return s.substring(fromAndTo[0], fromAndTo[1]);		
	}
	
	/**
	 * Computes common prefix of suffixes of two strings in O(N1 + N2 + M),
	 * N1 & N2 are length of strings and M is the number of LCP queries.
	 *
	 * <pre>
	 * LongestCommonPrefix lcp = new LongestCommonPrefix("acgat", "cgt");
	 * int[][] lcpQueries = new int[][]{
	 *     {1, 0}, // get common prefix length of "cgat" and "cgt"
	 *     {2, 1}, // get common prefix length of "gat" and "gt"
	 * };
	 * LongestCommonPrefix.Result lcpResult = lcp.preprocess(
	 *     Arrays.asList(lcpQueries));
	 * lcpResult.get(lcpQueries[0]) // returns 2 (common prefix is "cg")
	 * lcpResult.get(lcpQueries[1]) // returns 1 (common prefix is "c")
	 * </pre>
	 * @author bsaberid
	 *
	 */
	static class LongestCommonPrefix {
		/**
		 * Root of suffix tree
		 */
		private final Node root;
		/**
		 * For each string (s1 & s2), maintain a suffix table which store
		 * leaves in the reverse order of suffix length,
		 * e.g. get(0)[3] is the leaf node representing the 3rd suffix of s1,
		 * get(1)[0] is the leaf node representing the whole string s2.
		 */
		private final List<Node[]> suffixTable;
		/**
		 * Map which stores suffix length for each node in suffix tree
		 */
		private final Map<Node, Integer> depthMap;
		
		public LongestCommonPrefix(String s1, String s2) {
			StringBuilder sb = new StringBuilder(
					s1.length() + s2.length() + 2);
			sb.append(s1).append('\1').append(s2).append('\2');
			int[] endMarkerPos = new int[]{0, s1.length() + 1, sb.length()};
			
			root = SuffixTree.buildSuffixTree(sb);
			depthMap = new HashMap<Node, Integer>();
			suffixTable = new ArrayList<Node[]>(2);
			for (int i = 0; i < 2; ++i) {
				suffixTable.add(new Node[endMarkerPos[i + 1] -
				                   endMarkerPos[i]]);
			}
			dfs(root, endMarkerPos[1], suffixTable, depthMap);
		}
		
		private static final class NodeAndInt {
			final Node node; final int i;
			NodeAndInt(Node n, int i) {
				this.node = n;
				this.i = i;
			}
		}

		// we keep it static and have output parameters so we can test it easily
		private static void dfs(Node root,
				int firstEndMarkerPos,
				List<Node[]> suffixTable,
				Map<Node, Integer> depthMap) {
			ArrayDeque<NodeAndInt> stack = new ArrayDeque<NodeAndInt>();
			stack.add(new NodeAndInt(root, 0));
			while (!stack.isEmpty()) {
				NodeAndInt ni = stack.pollFirst();
				Node node = ni.node; int depth = ni.i;

				// compute the suffix length for each node in depthMap
				// (the SuffixTree.Node depth is node's depth in the tree) 
				int end = (node.begin > firstEndMarkerPos ? node.end 
						: Math.min(node.end, firstEndMarkerPos));
				depth += end - node.begin;
				depthMap.put(node, depth);

				boolean nonLeaf = false;
				for (Node child : node.children) {
					if (child != null) {
						nonLeaf = true;
						stack.addFirst(new NodeAndInt(child, depth));
					}
				}

				if (!nonLeaf && depth != 0) { // if leaf (and non-root if empty)
					int whichString = (node.begin <= firstEndMarkerPos ? 0 : 1);
					// we store leaves in the reverse order of depth
					int n = suffixTable.get(whichString).length;
					// n - (depth - 1) - 1 == n - depth
					suffixTable.get(whichString)[n - depth] = node;
				}
			}
		}
		
		// this method does not modify the state of LongestCommonPrefix class,
		// so can be used in parallel from different threads (each will have
		// its own result)
		public Result preprocess(Collection<int[]> suffixPosPairs) {
			Map<Node, NodeData> dataMap = new HashMap<Node, NodeData>(
					depthMap.size());
			for (Node node : depthMap.keySet())
				dataMap.put(node, new NodeData());
			Node[] firstSuffixTable = suffixTable.get(0);
			Node[] secondSuffixTable = suffixTable.get(1);
			for (int[] suffixPosPair : suffixPosPairs) {
				Node firstNode = firstSuffixTable[suffixPosPair[0]];
				Node secondNode = secondSuffixTable[suffixPosPair[1]];
				dataMap.get(firstNode).secondNodes.add(secondNode);
				dataMap.get(secondNode).secondNodes.add(firstNode);
			}
			
			Map<LCAQuery, Node> lcaMap = new HashMap<LCAQuery, Node>(
					suffixPosPairs.size()); // preallocate the right size
			tarjanOfflineLCA(root, dataMap, lcaMap);
			return new Result(lcaMap);
		}
		
		/**
		 * Proxy which returns preprocessed results in O(1) time.
		 */
		public class Result {
			private final Map<LCAQuery, Node> lcaMap;
			Result(Map<LCAQuery, Node> lcaMap) { this.lcaMap = lcaMap; }
			public int get(int[] suffixPos) {
				// can be extended to work with many nodes 					
				Node[] nodes = new Node[2];
				for (int i = 0; i < nodes.length; ++i)
				  nodes[i] = suffixTable.get(i)[suffixPos[i]];
				Node lca = lcaMap.get(new LCAQuery(nodes[0], nodes[1]));
				return depthMap.get(lca);
			}
		}
		
		private static class NodeData {
			boolean black;
			DisjointSetNode ds;
			List<Node> secondNodes = new ArrayList<Node>();
		}
		
		private static final class LCAQuery {
			final Node first, second; // Node[2] is more memory (has length)
			public LCAQuery(Node first, Node second) {
				this.first = first;
				this.second = second;
			}
			@Override
			public int hashCode() { // must be same if we swap first & second
				return first.hashCode() ^ second.hashCode();
			}
			@Override
			public boolean equals(Object obj) {
				LCAQuery other = (LCAQuery) obj;
				return (first == other.first && second == other.second ||
						first == other.second && second == other.first);
			}
		}
		
		private static void tarjanOfflineLCA(Node u,
				Map<Node, NodeData> dataMap,
				Map<LCAQuery, Node> lcaMap) {
			// This recursive version throws StackOverflowError for deep trees:
//			DisjointSetNode uSet = new DisjointSetNode(u);
//			TarjanOfflineData uData = dataMap.get(u);
//			uData.ds = uSet;
//			for (Node v : u.children) {
//				if (v == null) continue;
//				tarjanOfflineLCA(v, dataMap, lcaMap);
//				uSet.union(dataMap.get(v).ds);
//				uSet.find().ancestor = u;
//			}
//			uData.black = true;
//			for (Node v : uData.secondNodes) {
//				TarjanOfflineData vData = dataMap.get(v);
//				if (vData.black) {
//					Node lca = DisjointSetNode.find(vData.ds).ancestor;
//					lcaMap.put(new NodePair(u, v), lca);
//				}
//			}

			// So, we implement an iterative version:
			// - we keep track of child index which was visited last
			// - that way we can return to u after we visit each child
			ArrayDeque<NodeAndInt> stack = new ArrayDeque<NodeAndInt>();
			stack.push(new NodeAndInt(u, -1)); // -1 means no children visited
			
			while (!stack.isEmpty()) {
				NodeAndInt state = stack.pollFirst();
				u = state.node; int vIndex = state.i;
				NodeData uData = dataMap.get(u);
				if (vIndex < 0) {
					uData.ds = new DisjointSetNode(u);
				} else {
					DisjointSetNode uSet = uData.ds;
					DisjointSetNode vSet = dataMap.get(u.children[vIndex]).ds;
					uSet.union(vSet);
					uSet.find().ancestor = u;
				}
				
				boolean childrenVisited = true;
				while (++vIndex < u.children.length) {
					Node v = u.children[vIndex];
					if (v == null) continue;
					stack.addFirst(new NodeAndInt(u, vIndex)); // back to u
					stack.addFirst(new NodeAndInt(v, -1)); // recurse on v
					childrenVisited = false;
					break;
				}
				if (!childrenVisited)
					continue;
				
				uData.black = true;
				for (Node v : uData.secondNodes) {
					NodeData vData = dataMap.get(v);
					if (vData.black) {
						Node lca = vData.ds.find().ancestor;
						lcaMap.put(new LCAQuery(u, v), lca);
					}
				}
			}
		}
		
		static class DisjointSetNode {
			DisjointSetNode parent;
			int rank;
			Node ancestor;
			
			public DisjointSetNode(Node ancestor) {
				parent = this;
				rank = 0;
				this.ancestor = ancestor;
			}
			
			void union(DisjointSetNode other) {
				DisjointSetNode xRoot = find();
				DisjointSetNode yRoot = other.find();
				if (xRoot.rank < yRoot.rank)
					xRoot.parent = yRoot;
				else {
					yRoot.parent = xRoot;
					if (yRoot.rank == xRoot.rank)
						xRoot.rank += 1;
				}
			}
			
			DisjointSetNode find() {
				if (parent == this)
					return this;
				parent = parent.find(); // relink to root
				return parent; // return root
			}
		}
	}
}
