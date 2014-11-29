package longestPalindrome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import longestPalindrome.SuffixTree.Node;
import static org.junit.Assert.*;

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
	static void assertLongestPali(String s, String explongestPali) {
		//assert false, true : prevent the repetition
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
		
		public LongestCommonPrefix(String s1, String s2) {
			String[] s = new String[]{s1, s2};
			StringBuilder sb = new StringBuilder(
					s[0].length() + s[1].length() + 2);
			sb.append(s[0]).append('\1').append(s[1]).append('\2');
			int[] endMarkerPos = new int[]{0, s[0].length() + 1, sb.length()};
			
			root = SuffixTree.buildSuffixTree(sb);
			depthMap = new HashMap<Node, Integer>();
			suffixTable = new ArrayList<Node[]>(2);
			for (int i = 0; i < 2; ++i) {
				suffixTable.add(new Node[endMarkerPos[i + 1] -
				                   endMarkerPos[i]]);
			}
			dfs(root, 0, 0, endMarkerPos[1], suffixTable, depthMap);
			
			for (Node[] nodes : suffixTable) {
				System.err.println("New String");
				for (Node n : nodes)
					System.err.println("Node: " + n);
			}
		}
		
		private static void dfs(Node node, int index, int depth,
				int firstEndMarkerPos, 
				List<Node[]> suffixTable,
				Map<Node, Integer> depthMap) {
			if (node == null)
				return;
			
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
						firstEndMarkerPos, suffixTable, depthMap);
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
		
		public Node[] getLCAs(List<int[]> suffixPosPairs) {
			// TODO
		}
		
		private static class TarjanOfflineData {
			DisjointSetNode ds;
			boolean black;
			List<Node> secondNodes;
			List<Node> lcas;
		}
		
		static class NodePair {
			Node first, second;
			public NodePair(Node first, Node second) {
				this.first = first;
				this.second = second;
			}			
		}
		
		private DisjointSetNode tarjanOfflineLCA(Node u,
				Map<Node, TarjanOfflineData> dataMap,
				Map<NodePair, Node> lcaMap) {

			DisjointSetNode uSet = new DisjointSetNode(u);
			TarjanOfflineData uData = dataMap.get(u);
			uData.ds = uSet;
			for (Node v : u.children) {
				if (v == null) continue;
				DisjointSetNode vSet = tarjanOfflineLCA(v, dataMap);
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
	}
	
	static String findLongestPalindrome(String s) {
		LongestCommonPrefix lcp = new LongestCommonPrefix(
				s,
				new StringBuilder(s).reverse().toString());
		int n = s.length();
		int maxPaliLen = 1;
		int[] fromAndTo = new int[]{0, 1};
		
		List<int[]> lcaOddQuery = new ArrayList<int[]>(n);
		List<int[]> lcaEvenQuery = new ArrayList<int[]>(n);
		for (int i = 0; i < n; ++i) {
			lcaOddQuery.add(new int[]{i, n - i - 1}); // precompute for odd
			lcaEvenQuery.add(new int[]{i, n - i}); // precompute for even
		}
		Node[] lcaOddResult = lcp.getLCAs(lcaOddQuery); // O(n)
		Node[] lcaEvenResult = lcp.getLCAs(lcaEvenQuery); // O(n)
		
		for (int i = 0; i < n; ++i) {
			// check for odd-length palindrome centered at pos i in O(1)
			int k = lcp.depthMap.get(lcaOddResult[i]);
			int paliLen = 2 * k - 1;
			if (paliLen > maxPaliLen) {
				paliLen = maxPaliLen;
				fromAndTo[0] = i - k + 1;
				fromAndTo[1] = i + k;
			}
			// check for even-length palindrome centered at pos i in O(1)
			k = lcp.depthMap.get(lcaEvenResult[i]);
			paliLen = 2 * k;
			if (paliLen > maxPaliLen) {
				paliLen = maxPaliLen;
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
		assertLongestPali("abba", "abba");
		assertLongestPali("salaabcbasala", "abcba");
	}
	
	@Test
	public void testSpace() {
		assertLongestPali("sallaa abcba alsa", "abcba"); // check for the spaces
		assertLongestPali("madamis a doctor", "madam"); 
	}

	@Test	
	public void testSinglelong() {
		assertLongestPali("aaabmbaaadokhtar", "aaabmbaaa");      //beginning
		assertLongestPali("fghijaabmbaakmdafr", "aabmbaa");     //middle
		assertLongestPali("dokhtareaaaacbbbcaa", "aacbbbcaa");  //end

	}
	
	@Test
	public void testSingleShort() {
		assertLongestPali("salamdokhtarezibaabbasakamsdfchetori", "abba"); //beginning
		assertLongestPali("salamdokhtarezibasakamsdfchetoriabba", "abba"); //middle
		assertLongestPali("salamabbadokhtarezibasakamsdfchetori", "abba"); //end		
	}
	
	@Test
	public void testLongest() {
		// in case of two palis, return the longest one
		assertLongestPali("madamis a nun", "madam");
		assertLongestPali("salamdokhtarezibabcbasakamsdfabbachetori", "abcba"); 	
		assertLongestPali("salamabbadokhtarezibabcbasakamsdfchetori", "abcba");		
		assertLongestPali("the racecar has a race radar to check where is racecaring!", "racecar");
		assertLongestPali("radarof the the racecar checks where is he", "racecar");
	}
		
	@Test
	public void testOrder() {
		assertLongestPali("abscdxxxxxxxxxxxxxasdfsbbbbbbbbdfyyy", "xxxxxxxxxxxxx"); //beginning
		assertLongestPali("abcdfe3333adsfjsdf22674903454444444asdfjl;adsfjadsf11yy", "4444444"); //middle
		assertLongestPali("abyyysdfsd22222fjdsfxxxxxxxx", "xxxxxxxx");	//end		 
	}
	
	@Test 
	public void testOverlap() {
		assertLongestPali("acabacab", "acabaca");
		assertLongestPali("abaacaaba", "abaacaaba"); // the longest one is the string itself
		assertLongestPali("abaacaabaabaacaabaabaacaaba", "abaacaabaabaacaabaabaacaaba");
	}
}


