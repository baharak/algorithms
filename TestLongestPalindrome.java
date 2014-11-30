package longestPalindrome;

import static org.junit.Assert.*;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import longestPalindrome.SuffixTree.Node;

import org.junit.Test;

public class TestLongestPalindrome {

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
	
	public static String findLongestPalindrome(String s) {
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
	
	// TEST CASES
	
	@Test
	public void testLongestCommonPrefix() {
		LongestCommonPrefix lcp = new LongestCommonPrefix("acgat", "cgt");
		 int[][] lcpQueries = new int[][]{
		      {1, 0}, // get common prefix length of "cgat" and "cgt"
		      {2, 1}, // get common prefix length of "gat" and "gt"
		  };
		  LongestCommonPrefix.Result lcpResult = lcp.preprocess(
				  Arrays.asList(lcpQueries));
		  assertEquals(2, lcpResult.get(lcpQueries[0])); // "cg"
		  assertEquals(1, lcpResult.get(lcpQueries[1])) ;// "c"
	}
	
	private static void assertLongestPali(String s, String expLongestPali) {
		//assert false, true : prevent the repetition
		assertEquals(expLongestPali, findLongestPalindrome(s));
	}
	
	@Test
	public void testEmpty(){
		assertLongestPali("", "");
	}
	
	@Test
	public void testOddLength() {
		assertLongestPali("x", "x");
		assertLongestPali("aba", "aba");
		assertLongestPali("madamisadoctor", "madam");
		assertLongestPali("salaabcbasala", "abcba");
		assertLongestPali("salonisblas", "s");
		assertEquals(findLongestPalindrome("salonisblas").length(), 1);
	}
	
	@Test
	public void testEvenLength() {
		assertLongestPali("salamdokhtarezibaabbasakamsdfchetori", "baab");
		assertLongestPali("salamdokhtarezibasakamsdfchetoriabba", "abba");
		assertLongestPali("abba", "abba");
		assertLongestPali("moon", "oo");
		assertLongestPali("allumnimulla", "ll");
		assertLongestPali("22", "22");
	}

	@Test	
	public void testPosition() {
		assertLongestPali("aaabmbaaadokhtar", "aaabmbaaa"); 	//beginning
		assertLongestPali("fghijaabmbaakmdafr", "aabmbaa"); 	//middle
		assertLongestPali("dokhtareaaaacbbbcaa", "aacbbbcaa");	//end
	}
	
	@Test
	public void testLongest() {
		// in case of two palis, return the longest one
		assertLongestPali("madamisanun", "madam");
		assertLongestPali("salamdokhtarezibabcbasakamsdfabbachetori", "abcba");	
		assertLongestPali("salamabbadokhtarezibabcbasakamsdfchetori", "abcba");		
		assertLongestPali("theracecarhasaraceradartocheckwhereisracecaring",
				"racecar");
		assertLongestPali("radarofthetheracecarcheckswhereishe", "racecar");
	}
		
	@Test
	public void testOrder() {
		assertLongestPali("abscdxxxxxxxxxxxxxasdfsbbbbbbbbdfyyy",
				"xxxxxxxxxxxxx"); //beginning
		assertLongestPali(
				"abcdfe3333adsfjsdf22674903454444444asdfjladsfjadsf11yy",
				"4444444"); //middle
		assertLongestPali("abyyysdfsd22222fjdsfxxxxxxxx", "xxxxxxxx"); //end		 
	}
	
	@Test 
	public void testOverlap() {
		assertLongestPali("acabacab", "acabaca");
		assertLongestPali("m4321m1234", "4321m1234");
		assertLongestPali("abbacabbaaabb", "abbacabba");
		assertLongestPali("bbaaabbacabba", "abbacabba");
	}
	
	private static StringBuilder rep(char c, int n) {
		StringBuilder sb = new StringBuilder(n);
		for (int i = 0; i < n; ++i) 
			sb.append(c);
		return sb;
	}
	
	private static StringBuilder concat(CharSequence... charSeqs) {
		int len = 0;
		for (CharSequence charSeq : charSeqs)
			len += charSeq.length();
		StringBuilder sb = new StringBuilder(len);
		for (CharSequence charSeq : charSeqs)
			sb.append(charSeq);
		return sb;
	}
	
	private static StringBuilder cycleLetters(int n) {
		StringBuilder sb = new StringBuilder(n);
		char letter = 'a';
		for (int i = 0; i < n; ++i) {
			sb.append(letter);
			if (++letter > 'z') letter = 'a';
		}
		return sb;
	}
	
	private static StringBuilder pali(CharSequence charSeq) {
		StringBuilder sb = new StringBuilder(charSeq.length() * 2);
		sb.append(charSeq).reverse().append(charSeq);
		return sb;
	}
	
	// Tests in which there is a long palindrom of repeated character;
	// this creates a deep suffix tree for which we need non-recursive algorithm
	// (the depth of the suffix tree is O(n))
	
	@Test
	public void testARepeated100k() {
		final int n = 100000;
		StringBuilder sb = rep('a', n);
		String longestPali = sb.toString();
		assertLongestPali(sb.toString(), longestPali);
	}

	@Test
	public void testARepeated33kOf100k() {
		final int n = 33333;
		StringBuilder sb = rep('a', n + 1);
		String longestPali = sb.toString();
		assertLongestPali(concat(rep('c', n), sb, rep('t', n)).toString(),
				longestPali); // ??? aaa ???
		assertLongestPali(sb.append(rep('b', n)).append(rep('c', n)).toString(),
				longestPali); // aaa ??? ???
		assertLongestPali(sb.reverse().toString(), longestPali); // ??? ??? aaa
	}
	
	// Tests in which the longest palindrome is in the middle of the
	// first or second half of the long string (so that both algorithm that
	// check from the middle or from begin/end is slow):
	
	@Test
	public void testLongPaliAt2ndTo3rdSixthOf100k() {
		// ??? a-z...z-a ??? ??? ??? ???
		final int n = 100000;
		final int nOver6 = n / 6;
		StringBuilder longestPali = pali(cycleLetters(nOver6 / 2 + 1));
		StringBuilder sb = concat(
				rep('1', nOver6), longestPali, rep('3', nOver6),
				rep('4', nOver6), rep('5', nOver6), rep('6', nOver6));
		assertLongestPali(sb.toString(), longestPali.toString());
	}
	
	@Test
	public void testShortPaliAt3FourthsOf100k() {
		final int n = 100000;
		final int nOver4 = n / 4;
		StringBuilder longestPali = pali(rep('x', 7));
		StringBuilder sb = concat( // ??? ??? ??? xxxxxxx ???
				cycleLetters(3 * nOver4), longestPali, cycleLetters(nOver4));
		assertTrue(sb.length() > n);
		assertLongestPali(sb.toString(), longestPali.toString());
	}
}
