package longestPalindrome;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import longestPalindrome.SuffixTree.Node;

import org.junit.Ignore;
import org.junit.Test;

public class TestLongestPalindrome {
	
	static String Palindrome(String s) {
		String longestPalindrome = null;
		return longestPalindrome;
	}
	
	
	/**
	 * 
	 * @param s
	 * @param explongestPali
	 */
	static void assertLongestPali(String s, String expLongestPali) {
		//assert false, true : prevent the repetition
		assertEquals(expLongestPali, findLongestPalindrome(s));
	}
	
	static void generateString(String pattern) {
		
	}
	
	@Deprecated
	static void printNode(SuffixTree.Node node) {
		System.err.println("[begin=" + node.begin + ", end=" + node.end + "]");
	}
	
	//static int leafCnt = 0;
	
	static class LongestCommonPrefix {
		final Node root;
		/**
		 * For each string (s1 & s2), maintain a suffix table which store
		 * leaves in the reverse order of suffix length,
		 * e.g. get(0)[3] is the leaf node representing the 3rd suffix of s1,
		 * get(1)[0] is the leaf node representing the whole string s2.
		 */
		final List<Node[]> suffixTable;
		
		final Map<Node, Integer> depthMap;
		
		final Map<Node, TarjanOfflineData> dataMap;
		
		public LongestCommonPrefix(String s1, String s2) {
			String[] s = new String[]{s1, s2};
			StringBuilder sb = new StringBuilder(
					s[0].length() + s[1].length() + 2);
			sb.append(s[0]).append('\1').append(s[1]).append('\2');
			int[] endMarkerPos = new int[]{0, s[0].length() + 1, sb.length()};
			
			root = SuffixTree.buildSuffixTree(sb);
			depthMap = new HashMap<Node, Integer>();
			dataMap = new HashMap<Node, TarjanOfflineData>();
			suffixTable = new ArrayList<Node[]>(2);
			for (int i = 0; i < 2; ++i) {
				suffixTable.add(new Node[endMarkerPos[i + 1] -
				                   endMarkerPos[i]]);
			}
			dfs(root, 0, 0, endMarkerPos[1], suffixTable, depthMap, dataMap);
		}
		
		// we keep it static and have output parameters so we can test it easily
		private static void dfs(Node node, int index, int depth,
				int firstEndMarkerPos, 
				List<Node[]> suffixTable,
				Map<Node, Integer> depthMap,
				Map<Node, TarjanOfflineData> dataMap) {
			if (node == null)
				return;
			
			TarjanOfflineData data = new TarjanOfflineData();
			data.secondNodes = new ArrayList<Node>();
			dataMap.put(node, data);
			
			// compute the depth based on the suffix length
			int end = (node.begin > firstEndMarkerPos ? node.end 
					: Math.min(node.end, firstEndMarkerPos));
			depth += end - node.begin;
			depthMap.put(node, depth);
			
			boolean nonLeaf = false;
			int childIndex = 0;
			for (Node child : node.children) {
				nonLeaf |= (child != null);
				dfs(child, childIndex++, depth,
						firstEndMarkerPos, suffixTable, depthMap, dataMap);
			}

			if (!nonLeaf && depth != 0) { // if leaf (and non-root if empty)
				int whichString = (node.begin <= firstEndMarkerPos ? 0 : 1);
				// we store leaves in the reverse order of depth
				int n = suffixTable.get(whichString).length;
				// n - (depth - 1) - 1 == n - depth
				suffixTable.get(whichString)[n - depth] = node;
				//System.err.println("LEAF #" + leafCnt++ + ": " + node +
				//        " depth=" + depth + "[" + node.begin + ", " + end);
			}
		}
		
		public Map<NodePair, Node> getLCAs(List<int[]> suffixPosPairs) {
			Node[] firstSuffixTable = suffixTable.get(0);
			Node[] secondSuffixTable = suffixTable.get(1);
			for (int[] suffixPosPair : suffixPosPairs) {
				Node firstNode = firstSuffixTable[suffixPosPair[0]];
				Node secondNode = secondSuffixTable[suffixPosPair[1]];
				dataMap.get(firstNode).secondNodes.add(secondNode);
				dataMap.get(secondNode).secondNodes.add(firstNode);
			}
			
			Map<NodePair, Node> lcaMap = new HashMap<NodePair, Node>();
			tarjanOfflineLCA(root, dataMap, lcaMap);
			return lcaMap;
		}
		
		private static class TarjanOfflineData {
			DisjointSetNode ds;
			boolean black;
			List<Node> secondNodes;
		}
		
		// Node[2] is more memory (has length)
		public static final class NodePair {
			final Node first, second;
			public NodePair(Node first, Node second) {
				this.first = first;
				this.second = second;
			}
			@Override
			public int hashCode() {
				return first.hashCode() ^ second.hashCode();
			}
			@Override
			public boolean equals(Object obj) {
				NodePair other = (NodePair) obj;
				return (first == other.first && second == other.second ||
						first == other.second && second == other.first);
			}
		}
		
		private static DisjointSetNode tarjanOfflineLCA(Node u,
				Map<Node, TarjanOfflineData> dataMap,
				Map<NodePair, Node> lcaMap) {

			DisjointSetNode uSet = new DisjointSetNode(u);
			TarjanOfflineData uData = dataMap.get(u);
			uData.ds = uSet;
			for (Node v : u.children) {
				if (v == null) continue;
				DisjointSetNode vSet = tarjanOfflineLCA(v, dataMap, lcaMap);
				DisjointSetNode.union(uSet, vSet);
				DisjointSetNode.find(uSet).ancestor = u;
			}
			uData.black = true;
			
			for (Node v : uData.secondNodes) {
				TarjanOfflineData vData = dataMap.get(v);
				if (vData.black) {
					Node lca = DisjointSetNode.find(vData.ds).ancestor;
					lcaMap.put(new NodePair(u, v), lca);
				}
			}
			
			return uSet;
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
			
			static void union(DisjointSetNode x, DisjointSetNode y) {
				DisjointSetNode xRoot = find(x);
				DisjointSetNode yRoot = find(y);
				if (xRoot.rank < yRoot.rank)
					xRoot.parent = yRoot;
				else {
					yRoot.parent = xRoot;
					if (yRoot.rank == xRoot.rank)
						xRoot.rank += 1;
				}
					
			}
			
			static DisjointSetNode find(DisjointSetNode x) {
				if (x.parent == x)
					return x;
				x.parent = find(x.parent); // relink to root
				return x.parent; // return root
			}
			
		}
	
		public int getLCADepth(int[] suffixPos,
				// can be extended to work with many nodes 
				Map<NodePair, Node> lcaMap) {
			Node[] nodes = new Node[2];
			for (int i = 0; i < nodes.length; ++i)
			  nodes[i] = suffixTable.get(i)[suffixPos[i]];
			Node lca = lcaMap.get(
					new LongestCommonPrefix.NodePair(nodes[0], nodes[1]));
//			System.err.println("QRY(" + nodes[0] + "," + nodes[1] + "): "
//					+ lcaMap.size());
			return depthMap.get(lca);
		}
	}
	
	static String findLongestPalindrome(String s) {
		int n = s.length();
		int maxPaliLen = 0;
		int[] fromAndTo = new int[]{0, 0};
		List<int[]> lcaOddQuery = new ArrayList<int[]>(n);
		List<int[]> lcaEvenQuery = new ArrayList<int[]>(n);
		for (int i = 0; i < n; ++i) {
			lcaOddQuery.add(new int[]{i, n - i - 1}); // precompute for odd
			lcaEvenQuery.add(new int[]{i, n - i}); // precompute for even
		}
		
		LongestCommonPrefix lcp = new LongestCommonPrefix(
				s,
				new StringBuilder(s).reverse().toString());
		// TODO: lcp.preprocess(lcaQueries)
		Map<LongestCommonPrefix.NodePair, Node> lcaOdd =
				lcp.getLCAs(lcaOddQuery); // O(n)
		Map<LongestCommonPrefix.NodePair, Node> lcaEven =
				lcp.getLCAs(lcaEvenQuery); // O(n)
		for (int i = 0; i < n; ++i) {
			// check for odd-length palindrome centered at pos i in O(1)
			int k = lcp.getLCADepth(lcaOddQuery.get(i), lcaOdd);
			int paliLen = 2 * k - 1;
//			System.err.println("LCP(" 
//					+ Arrays.toString(lcaOddQuery.get(i)) + ") = " + k);
			if (paliLen > maxPaliLen) {
				maxPaliLen = paliLen;
				fromAndTo[0] = i - k + 1;
				fromAndTo[1] = i + k;
			}
			// check for even-length palindrome centered at pos i in O(1)
			k = lcp.getLCADepth(lcaEvenQuery.get(i), lcaEven);
			paliLen = 2 * k;
			if (paliLen > maxPaliLen) {
				maxPaliLen = paliLen;
				fromAndTo[0] = i - k;
				fromAndTo[1] = i + k;
			}
		}
		return s.substring(fromAndTo[0], fromAndTo[1]);		
	}
	
	@Ignore
	public void testSuffixTreeHugeString() {
		StringBuilder sb = new StringBuilder("acacag");
		for (int i = 0; i < 17; ++i)
			sb.append(sb.toString());
		assertTrue(sb.length() > 600000);
		String s1 = sb.toString();
		String s2 = sb.reverse().toString();
		SuffixTree.buildSuffixTree(s1 + '\1' + s2 + '\2');
	}
	
	@Test
	public void testSuffixTreeLeafsByDepth() {
		String s1 = "acgat", s2 = "cgt";
		LongestCommonPrefix lcp = new LongestCommonPrefix(s1, s2);
	}
	
	@Test
	public void testNull(){
		assertLongestPali("", "");
	}
	
	@Test
	public void testEasy() {
		assertLongestPali("aba", "aba");
		//assertLongestPali("abba", "abba");
		assertLongestPali("salaabcbasala", "abcba");
	}
	
	@Test
	public void testSpace() {
		assertLongestPali("sallaaabcbaalsa", "aabcbaa"); // check for the spaces
		assertLongestPali("madamisadoctor", "madam");
	}

	@Test	
	public void testSinglelong() {
		assertLongestPali("aaabmbaaadokhtar", "aaabmbaaa");      //beginning
		assertLongestPali("fghijaabmbaakmdafr", "aabmbaa");     //middle
		assertLongestPali("dokhtareaaaacbbbcaa", "aacbbbcaa");  //end

	}
	
	@Test
	public void testSingleShort() {
		assertLongestPali("salamdokhtarezibaabbasakamsdfchetori", "baab"); //beginning
		assertLongestPali("salamdokhtarezibasakamsdfchetoriabba", "abba"); //middle
		assertLongestPali("salamabbadokhtarezibasakamsdfchetori", "abba"); //end		
	}
	
	@Test
	public void testLongest() {
		// in case of two palis, return the longest one
		assertLongestPali("madamisanun", "madam");
		assertLongestPali("salamdokhtarezibabcbasakamsdfabbachetori", "abcba");	
		assertLongestPali("salamabbadokhtarezibabcbasakamsdfchetori", "abcba");		
		assertLongestPali("theracecarhasaraceradartocheckwhereisracecaring", "racecar");
		assertLongestPali("radarofthetheracecarcheckswhereishe", "racecar");
	}
		
	@Test
	public void testOrder() {
		assertLongestPali("abscdxxxxxxxxxxxxxasdfsbbbbbbbbdfyyy", "xxxxxxxxxxxxx"); //beginning
		assertLongestPali("abcdfe3333adsfjsdf22674903454444444asdfjladsfjadsf11yy", "4444444"); //middle
		assertLongestPali("abyyysdfsd22222fjdsfxxxxxxxx", "xxxxxxxx");	//end		 
	}
	
	@Test 
	public void testOverlap() {
		assertLongestPali("acabacab", "acabaca");
		assertLongestPali("abaacaaba", "abaacaaba"); // the longest one is the string itself
		assertLongestPali("abaacaabaabaacaabaabaacaaba", "abaacaabaabaacaabaabaacaaba");
	}
}


