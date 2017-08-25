﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class DatabaseReference
    {
        #region Public Constructors

        /// <summary>
        /// DatabaseReference constructor, which takes the type and ID strings of the reference, and a list of properties. Each property contains the "type" and "value" of the property as Item1 and Item2 of the Tuple.
        /// </summary>
        /// <param name="type"></param>
        /// <param name="id"></param>
        /// <param name="properties"></param>
        public DatabaseReference(string type, string id, IEnumerable<Tuple<string, string>> properties)
        {
            Type = type;
            Id = id;
            Properties = properties;
        }

        #endregion Public Constructors

        #region Public Properties

        /// <summary>
        /// dbRef type, e.g. "GO" for GO terms
        /// </summary>
        public string Type { get; }

        /// <summary>
        /// dbRef ID string
        /// </summary>
        public string Id { get; }

        /// <summary>
        /// Each database reference contains a list of properties. Item1 of this Tuple is the "type", and Item2 is the "value" of the property.
        /// </summary>
        public IEnumerable<Tuple<string, string>> Properties { get; }

        #endregion Public Properties

        #region Public Methods

        public override bool Equals(object obj)
        {
            DatabaseReference d = obj as DatabaseReference;
            return obj != null
                && d.Type == Type
                && d.Id == Id
                && d.Properties.OrderBy(x => x).SequenceEqual(Properties.OrderBy(x => x));
        }

        public override int GetHashCode()
        {
            int hash = Type.GetHashCode() ^ Id.GetHashCode();
            foreach (Tuple<string, string> property in Properties)
            {
                hash = hash ^ property.GetHashCode();
            }
            return hash;
        }

        #endregion Public Methods
    }
}