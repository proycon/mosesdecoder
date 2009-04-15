#pragma once

#include <vector>

namespace Moses
{
class Word;

class CellCollection
{
public:
	virtual const std::vector<Word> &GetHeadwords(const Moses::WordsRange &coverage) const = 0;

};

}

