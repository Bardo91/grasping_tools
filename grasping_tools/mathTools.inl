///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mathTools.h"

namespace gpisGrasping {
	template<typename _vectorType>
	inline std::vector<_vectorType> generateCombinationsVaryingSize(const _vectorType &_vec) {
		std::vector<_vectorType> combinations;

		std::function<void(_vectorType, _vectorType)> recursiveAdd = [&](_vectorType _set, _vectorType _head) {
			for (unsigned i = 0; i < _set.size(); i++) {
				_vectorType toAdd = _head;
				toAdd.insert_rows(toAdd.size(), _set.row(i));
				combinations.push_back(toAdd);
			}
		
			for (unsigned i = 0; i < _set.size() - 1; i++) {
				_vectorType newHead = _head;
				newHead.insert_rows(newHead.size(), _set.row(i));
				_vectorType newSet = _set;
				newSet.shed_rows(0, i);
				recursiveAdd(newSet, newHead);
			}
		
		};
		
		recursiveAdd(_vec, {});

		return combinations;
	}

	//--------------------------------------------------------------------------------------------------------------------
	template<typename _vectorType>
	std::vector<_vectorType> generateCombinationsWithoutRepetitions(const _vectorType & _vec, unsigned k) {
		std::vector<_vectorType> combinations;
		
		std::function<void(_vectorType, _vectorType)> recursiveAdd = [&](_vectorType _head, _vectorType _set) {
			if (_set.size() == 0) {
				// Reject.
			}
			else if (_set.size() > 0 && _head.size() == k - 1) {	// Add combination.
				for (unsigned i = 0; i < _set.size(); i++) {
					_vectorType newHead = _head;
					newHead.insert_rows(_head.size(), _set.row(i));
					combinations.push_back(newHead);
				}
			}
			else {	// More recursivity.
				for (unsigned i = 0; i < _set.size() - 1 ; i++) {
					_vectorType newHead = _head;
					_vectorType newSet = _set.rows(i + 1, _set.size()-1);
					newHead.insert_rows(_head.size(), _set.row(i));
					recursiveAdd(newHead, newSet);
				}
			}


		};

		recursiveAdd({}, _vec);


		return combinations;
	}
}