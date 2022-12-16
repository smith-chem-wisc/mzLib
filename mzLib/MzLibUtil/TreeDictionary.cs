using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace MzLibUtil
{
    // Implemented as a basic binary search tree, where nodes are ordered according to 
    // a key double (generally, an m/z value) and a value (generally, the corresponding intensity value)
    // This affords efficient search, deletion, and insertion ( O(logn) )
    // Keeping the "key" as a double allows operations such as FindNearest
    public class SpectrumTree : IEnumerable
    {
        private Node _root;
        public Node Root => _root;

        public SpectrumTree()
        {
            _root = null;
        }

        public void Insert(Node newNode)
        {
            if (Root == null)
            {
                _root = newNode;
                return;
            }

            Node y = null;
            Node x = Root;
            while (x != null)
            {
                y = x;
                x = newNode.Key < x.Key ? x.leftChild : x.rightChild;
            }
            newNode.parent = y;

            if (newNode.Key < y.Key)
            {
                y.leftChild = newNode;
            }
            else
            {
                y.rightChild = newNode;
            }
        }

        public void Delete(Node node)
        {
            if (node.parent == null)
            {
                _root = null;
                return;
            }

            if (node.parent.leftChild == node)
            {
                // Leaf deletion
                if (node.leftChild == null & node.rightChild == null)
                {
                    node.parent.leftChild = null;
                    return;
                }
                // One child deletion
                if (node.leftChild == null & node.rightChild != null)
                {
                    node.parent.leftChild = node.rightChild;
                    node.rightChild.parent = node.parent;
                    return;
                }
                if (node.leftChild != null & node.rightChild == null)
                {
                    node.parent.leftChild = node.leftChild;
                    node.leftChild.parent = node.parent;
                    return;
                }
                // If the node has two children, swap values with the succesor, then delete the successor
                Node inOrderSuccesor = GetInOrderSuccessor(node);
                node.SwapKeyValuePair(inOrderSuccesor.Key, inOrderSuccesor.Value);
                Delete(inOrderSuccesor);
                
            } else if (node.parent.rightChild == node)
            {
                if (node.leftChild == null & node.rightChild == null)
                {
                    node.parent.rightChild = null;
                    return;
                }
                if (node.leftChild == null & node.rightChild != null)
                {
                    node.parent.rightChild = node.rightChild;
                    node.rightChild.parent = node.parent;
                    return;
                }
                if (node.leftChild != null & node.rightChild == null)
                {
                    node.parent.rightChild = node.leftChild;
                    node.leftChild.parent = node.parent;
                    return;
                }
                Node inOrderSuccesor = GetInOrderSuccessor(node);
                node.SwapKeyValuePair(inOrderSuccesor.Key, inOrderSuccesor.Value);
                Delete(inOrderSuccesor);
            }
        }

        internal static Node GetInOrderSuccessor(Node predeccesor)
        {
            if (predeccesor.rightChild != null)
            {
                return GetMinNode(predeccesor.rightChild);
            } 

            Node parent = predeccesor.parent;
            while (parent != null && predeccesor != parent.leftChild)
            {
                predeccesor = parent;
                parent = parent.parent;
            }

            return parent;
        }

        internal static Node GetMinNode(Node parent)
        {
            Node current = parent;

            while (current.leftChild != null)
            {
                current = current.leftChild;
            }

            return current;
        }

        public IEnumerator GetEnumerator()
        {
            return new TreeEnumerator(_root);
        }

        /// <summary>
        /// This function assumes, but does not require, that the mzArray is ordered.
        /// For ordered arrays, creates a balanced binary tree
        /// Unordered arrays are fine, but then the resulting tree will not be balanced
        /// </summary>
        /// <param name="mzArray"></param>
        /// <param name="intensityArray"></param>
        public void BuildTreeFromSpectrum(double[] mzArray, double[] intensityArray)
        {
            BuildTreeHelper(mzArray, intensityArray, 0, mzArray.Length-1);
        }

        private void BuildTreeHelper(double[] mzArray, double[] intensityArray, int minIndex, int maxIndex)
        {
            // Base case - the array cannot be subdivided further
            if (minIndex == maxIndex)
            {
                Insert( new Node(mzArray[minIndex], intensityArray[minIndex]) );
                return;
            }
            if (maxIndex < minIndex) return; // Other base case, happens on right hand side of subsection

            int midIndex = minIndex + (maxIndex-minIndex)/2;
            Insert(new Node(mzArray[midIndex], intensityArray[midIndex]));
            BuildTreeHelper(mzArray, intensityArray, minIndex, midIndex - 1);
            BuildTreeHelper(mzArray, intensityArray, midIndex + 1, maxIndex);
        }
    }

    public class Node
    {
        // Theoretically the key could be any IComparable, but then FindNearest wouldn't be possible
        private double _key;
        public double Key => _key;
        private double _value;
        public double Value => _value;

        internal Node parent;
        internal Node leftChild;
        internal Node rightChild;
        internal double inOrderSuccesor;
        internal double inOrderPredecessor;

        public Node(double key, double value)
        {
            _key = key;
            _value = value;
        }

        internal void SwapKeyValuePair(double newKey, double newValue)
        {
            _key = newKey;
            _value = newValue;
        }

        public override string ToString()
        {
            return Key + "," + Value;
        }
    }

    public class TreeEnumerator : IEnumerator
    {
        private Node _currentNode = null;
        private Node _nextNode = null;
        private Node _rootNode;

        public TreeEnumerator(Node root)
        {
            _rootNode = root;
            _nextNode = SpectrumTree.GetMinNode(_rootNode);
        }

        public bool MoveNext()
        {
            if (_nextNode != null)
            {
                _currentNode = _nextNode;
                _nextNode = SpectrumTree.GetInOrderSuccessor(_currentNode);
                return true;
            }

            return false;
        }

        public void Reset()
        {
            _currentNode = null;
            _nextNode = SpectrumTree.GetMinNode(_rootNode);
        }

        object IEnumerator.Current => _currentNode;
        public Node Current => _currentNode;

    }
}
