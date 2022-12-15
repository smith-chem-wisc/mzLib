using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace MzLibUtil
{
    // Implemented as a basic 
    public class TreeDictionary
    {
        private Node root;

        internal class Node
        {
            // Theoretically the key could be any IComparable, but then FindNearest wouldn't be possible
            private double _key;
            internal double Key => _key;
            private Object _value;
            public Object Value => _value;

            internal Node parent;
            internal Node leftChild;
            internal Node rightChild;
            internal double inOrderSuccesor;
            internal double inOrderPredecessor;

            public Node(double key, Object value)
            {
                _key = key;
                _value = value;
            }

            internal void SwapKeyValuePair(double newKey, Object newValue)
            {
                _key = newKey;
                _value = newValue;
            }
        }

        private void Insert(Node newNode)
        {
            if (root == null)
            {
                root = newNode;
                root.color = Color.Black;
                return;
            }

            Node y = null;
            Node x = root;
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

        private void Delete(Node node)
        {
            if (node.parent == null)
            {
                root = null;
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

        private Node GetInOrderSuccesor(Node predeccesor)
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

        private Node GetMinNode(Node parent)
        {
            Node current = parent;

            if (current.leftChild != null)
            {
                current = current.leftChild;
            }

            return current;
        }

        public void BuildTreeFromSpectrum(double[] mzArray, double[] intensityArray)
        {
            int middleIndex = mzArray.Length / 2;

        }

        private void BuildTreeHelper(double[] mzArray, double[] intensityArray, int middleIndex)
        {

        }
    }
}
