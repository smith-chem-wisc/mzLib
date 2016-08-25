// Copyright 2016 Stefan Solntsev
//
// This file (PeriodicTableValidationResult.cs) is part of Chemistry Library.
//
// Chemistry Library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chemistry Library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Chemistry Library. If not, see <http://www.gnu.org/licenses/>.

namespace Chemistry
{
    public enum ValidationResult
    {
        PassedAbundanceValidation,
        PassedAverageMassValidation,
        FailedAbundanceValidation,
        FailedAverageMassValidation,
    }

    public class PeriodicTableValidationResult
    {
        public ValidationResult ThisValidationResult { get; private set; }
        public string Message { get; private set; }

        public PeriodicTableValidationResult(ValidationResult validationResult, string message)
        {
            ThisValidationResult = validationResult;
            Message = message;
        }
    }
}