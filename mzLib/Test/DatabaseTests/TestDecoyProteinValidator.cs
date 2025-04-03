using NUnit.Framework;
using UsefulProteomicsDatabases.DecoyGeneration;

namespace Test.DatabaseTests;

[TestFixture]
public class DecoySequenceValidatorTests
{
    [TestCase("AABBAA", 3, 3, true, TestName = "Full Palindrome-Even Sequence Length-Cutoff Equal")]
    [TestCase("AABBAA", 3, null, true, TestName = "Full Palindrome-Even Sequence Length-Cutoff Null")]
    [TestCase("AABBAA", 3, 4, false, TestName = "Full Palindrome-Even Sequence Length-Cutoff Greater")]
    [TestCase("AABBAA", 3, 2, true, TestName = "Full Palindrome-Even Sequence Length-Cutoff Lesser")]
    [TestCase("ABCBA", 3, 4, false, TestName = "Full Palindrome-Odd Sequence Length-Cutoff Greater")]
    [TestCase("ABCBA", 3, null, true, TestName = "Full Palindrome-Odd Sequence Length-Cutoff Null")]
    [TestCase("ABCBA", 3, 3, true, TestName = "Full Palindrome-Odd Sequence Length-Cutoff Equal")]
    [TestCase("ABCBA", 3, 2, true, TestName = "Full Palindrome-Odd Sequence Length-Cutoff Lesser")]
    [TestCase("ABCDEFCBA", 3, 2, true, TestName = "Partial Palindrome-Odd Sequence-Cutoff Lesser")]
    [TestCase("ABCDEFCBA", 3, null, false, TestName = "Partial Palindrome-Odd Sequence-Cutoff Null")]
    [TestCase("ABCDEFCBA", 3, 3, true, TestName = "Partial Palindrome-Odd Sequence-Cutoff Equal")]
    [TestCase("ABCDEFCBA", 3, 4, false, TestName = "Partial Palindrome-Odd Sequence-Cutoff Greater")]
    [TestCase("ABCDEGFCBA", 3, 2, true, TestName = "Partial Palindrome-Even Sequence-Cutoff Lesser")]
    [TestCase("ABCDEGFCBA", 3, null, false, TestName = "Partial Palindrome-Even Sequence-Cutoff Null")]
    [TestCase("ABCDEGFCBA", 3, 3, true, TestName = "Partial Palindrome-Even Sequence-Cutoff Equal")]
    [TestCase("ABCDEGFCBA", 3, 4, false, TestName = "Partial Palindrome-Even Sequence-Cutoff Greater")]
    [TestCase("ABCDEGF", 0, 2, false, TestName = "Non-Palindrome-Odd Sequence-Cutoff Lesser")]
    [TestCase("ABCDEGF", 0, null, false, TestName = "Non-Palindrome-Odd Sequence-Cutoff Null")]
    [TestCase("ABCDEGF", 0, 3, false, TestName = "Non-Palindrome-Odd Sequence-Cutoff Equal")]
    [TestCase("ABCDEGFD", 0, 4, false, TestName = "Non-Palindrome-Odd Sequence-Cutoff Greater")]
    [TestCase("ABCDEGFD", 0, 2, false, TestName = "Non-Palindrome-Even Sequence-Cutoff Lesser")]
    [TestCase("ABCDEGFD", 0, null, false, TestName = "Non-Palindrome-Even Sequence-Cutoff Null")]
    [TestCase("ABCDEGFD", 0, 3, false, TestName = "Non-Palindrome-Even Sequence-Cutoff Equal")]
    [TestCase("ABCDEGFD", 0, 4, false, TestName = "Non-Palindrome-Even Sequence-Cutoff Greater")]
    [TestCase("", 0, 1, false, TestName = "Empty Sequence")]
    [TestCase(null, 0, 1, false, TestName = "Null Sequence")]
    public void IsPalindromicTests(string input, int expectedDegree, int? degreeCutoff, bool expectedResult)
    {
        bool result = DecoySequenceValidator.IsPalindromic(input, out int degreeOfPalindromicity, degreeCutoff);
        Assert.That(result, Is.EqualTo(expectedResult));
        Assert.That(degreeOfPalindromicity, Is.EqualTo(expectedDegree));
    }
}