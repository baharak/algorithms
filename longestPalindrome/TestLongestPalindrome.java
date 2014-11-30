/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package longestPalindrome;

import static org.junit.Assert.*;

import java.util.Arrays;

import longestPalindrome.LongestPalindrome.LongestCommonPrefix;

import org.junit.Test;

public class TestLongestPalindrome {
	
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
		  assertEquals(1, lcpResult.get(lcpQueries[1])); // "c"
	}
	
	private static void assertLongestPali(String s, String expLongestPali) {
		//assert false, true : prevent the repetition
		assertEquals(expLongestPali, LongestPalindrome.findIn(s));
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
		assertEquals(LongestPalindrome.findIn("salonisblas").length(), 1);
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
