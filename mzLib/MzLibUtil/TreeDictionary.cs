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
                if (newNode.Key < x.Key)
                {
                    x = x.leftChild;
                }
                else
                {
                    x = x.rightChild;
                }
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
                // If the node has two children, swap values (maybe not the best way to do this?)
                Node inOrderSuccesor = GetInOrderSuccesor(node);
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
                Node inOrderSuccesor = GetInOrderSuccesor(node);
                node.SwapKeyValuePair(inOrderSuccesor.Key, inOrderSuccesor.Value);
                Delete(inOrderSuccesor);
            }
        }

        internal static Node GetInOrderSuccesor(Node predeccesor)
        {
            if (predeccesor.rightChild != null)
            {
                return GetMinNode(predeccesor.rightChild);
            } else if (predeccesor.parent.rightChild == predeccesor)
            {
                return null;
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

            if (current.leftChild != null)
            {
                current = current.leftChild;
            }

            return current;
        }

        public IEnumerator GetEnumerator()
        {
            return new TreeEnumerator(_root);
        }

        //public List<Node> GetNodeList()
        //{
        //    List<Node> allNodes = new();
        //    Node currentNode = GetMinNode(Root);
        //    while (currentNode != null)
        //    {
        //        allNodes.Add(currentNode);
        //        currentNode = GetInOrderSuccesor(currentNode);
        //    }
        //    return allNodes;
        //}

        public void BuildTreeFromSpectrum(double[] mzArray, double[] intensityArray)
        {
            int middleIndex = mzArray.Length / 2;
        }

        private void BuildTreeHelper(double[] mzArray, double[] intensityArray, int middleIndex)
        {

        }
    }

    public class Node
    {
        // Theoretically the key could be any IComparable, but then FindNearest wouldn't be possible
        private double _key;
        internal double Key => _key;
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
                _nextNode = SpectrumTree.GetInOrderSuccesor(_currentNode);
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
