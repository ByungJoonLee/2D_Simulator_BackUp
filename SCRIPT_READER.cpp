#include "SCRIPT_READER.h"

const char* SCRIPT_READER::REGEX_INCLUDE_STRING  = "\\s*#INCLUDE\\s*<\\s*([\\.\\w\\/]+)\\s*>\\s*";

const char* SCRIPT_READER::REGEX_BLOCKDATA_STRING       = "\\s*(\\w+)\\s*\\{\\s*(.*)\\s*\\}\\s*";
const char* SCRIPT_READER::REGEX_ASSIGN_BOOLEAN_STRING  = "\\s*(\\w+)\\s*=\\s*(true|false)\\s*";
const char* SCRIPT_READER::REGEX_ASSIGN_NUMBER_STRING   = "\\s*(\\w+)\\s*=\\s*(-*\\d*\\.*\\d+)\\s*";
const char* SCRIPT_READER::REGEX_ASSIGN_VECTOR2_STRING  = "\\s*(\\w+)\\s*=\\s*(\\(\\s*-*\\d*\\.*\\d+\\s*,\\s*-*\\d*\\.*\\d+\\s*\\))\\s*";

const char* SCRIPT_READER::REGEX_ASSIGN_STRING_STRING   = "\\s*(\\w+)\\s*=\\s*\"(\\s*[\\:\\\\.\\w\\/\\s\\-.*\\m\\d]+\\s*)\"\\s*";

const char* SCRIPT_READER::REGEX_ASSIGN_VECTOR2_VALUE_STRING  = "\\s*\\(\\s*(-*\\d*\\.*\\d+)\\s*,\\s*(-*\\d*\\.*\\d+)\\s*\\)\\s*";

const boost::regex SCRIPT_READER::REGEX_VALUE_NAME     = boost::regex("\\s*(\\w+)\\s*");
const boost::regex SCRIPT_READER::REGEX_INCLUDE        = boost::regex(REGEX_INCLUDE_STRING);
const boost::regex SCRIPT_READER::REGEX_BLOCKDATA      = boost::regex(REGEX_BLOCKDATA_STRING);
const boost::regex SCRIPT_READER::REGEX_ASSIGN_NUMBER  = boost::regex(REGEX_ASSIGN_NUMBER_STRING );
const boost::regex SCRIPT_READER::REGEX_ASSIGN_VECTOR2 = boost::regex(REGEX_ASSIGN_VECTOR2_STRING);
const boost::regex SCRIPT_READER::REGEX_ASSIGN_BOOLEAN = boost::regex(REGEX_ASSIGN_BOOLEAN_STRING);
const boost::regex SCRIPT_READER::REGEX_ASSIGN_STRING  = boost::regex(REGEX_ASSIGN_STRING_STRING);

const boost::regex SCRIPT_READER::REGEX_ASSIGN_VECTOR2_VALUE = boost::regex(REGEX_ASSIGN_VECTOR2_VALUE_STRING);

const SCRIPT_BLOCK* SCRIPT_READER::empty_block = new SCRIPT_BLOCK;
const SCRIPT_VALUE* SCRIPT_READER::empty_value = new SCRIPT_VALUE;

bool SCRIPT_READER::Initialize(const char* script)
{
	ReleaseScriptReader();

	FILE* file_script;
	if (!(file_script = fopen(script, "rt")))
	{
		return false;
	}

	script_file_names.push_back(string(script));

	// text_script is the contents in .sim file
	string text_script;
	text_script.clear();

	char char_script;
	while ((char_script = fgetc(file_script)) != EOF)
	{
		text_script += char_script;
	}
	
	fclose(file_script);

	// Include script and rewrite source
	string result;
	PreprocessIncludeSource(text_script.data(), &result);

	// Cut block at source, Call ParseRegexBlock per each block
	vector<string> trim;

	TrimSourceBlock(result.data(), &trim);

	for (int i = 0; i < (int)trim.size(); i++)
	{
		ParseRegexScriptBlock(trim[i].data());
	}

	return true;
}

void SCRIPT_READER::ReleaseScriptBlock(SCRIPT_BLOCK* block)
{
	if (!block)
	{
		return;
	}
	
	int i = 0;
	int size = (int)block->values.size();
	
	for (i = 0; i < size; i++)
	{
		delete block->values[i];
	}

	i = 0;
	size = (int)block->child_blocks.size();

	for (i = 0; i < size; i++)
	{
		ReleaseScriptBlock(block->child_blocks[i]);
	}

	block->child_blocks.clear();
	block->values.clear();
	
	delete block;
}

void SCRIPT_READER::ReleaseScriptReader()
{
	int i = 0;
	int s = (int)script_blocks.size();
	for (i = 0; i < s; i++)
	{
		ReleaseScriptBlock(script_blocks[i]);
	}

	script_blocks.clear();
	script_file_names.clear();
}

void SCRIPT_READER::ParseRegexScriptBlock(const char* script, SCRIPT_BLOCK* pBlock)
{
	SCRIPT_BLOCK* new_block = new SCRIPT_BLOCK();		// Creating new data block

	if (!pBlock)
	{
		script_blocks.push_back(new_block);
	}
	else
	{
		pBlock->child_blocks.push_back(new_block);
	}

	boost::cmatch matches;
	boost::regex_search(script, matches, REGEX_BLOCKDATA);
	new_block->block_name = string(matches[1].first, matches[1].second);	// Setting block name
	string data_string = string(matches[2].first, matches[2].second);

	// Cut child block and generate new Data Block
	vector<string> trim;
	TrimSourceBlock(data_string.data(), &trim);

	ParseRegexScriptValue(data_string.data(), SCRIPT_VALUE::TYPE_NUMBER, new_block);
	ParseRegexScriptValue(data_string.data(), SCRIPT_VALUE::TYPE_VECTOR2, new_block);
	ParseRegexScriptValue(data_string.data(), SCRIPT_VALUE::TYPE_BOOLEAN, new_block);
	ParseRegexScriptValue(data_string.data(), SCRIPT_VALUE::TYPE_STRING, new_block);

	// Call ParseRegexBlock() again, if child block is exist
	for (int i = 0; i < (int)trim.size(); i++)
	{
		ParseRegexScriptBlock(trim[i].data(), new_block);
	}
}

void SCRIPT_READER::ParseRegexScriptValue(const char* script, const boost::regex* ex, SCRIPT_BLOCK* pBlock)
{
	char* cp = (char*)script;			// Character point
	boost::cmatch matches;

	while (boost::regex_search(cp, matches, *ex))
	{
		SCRIPT_VALUE* data = new SCRIPT_VALUE();
		data->value_name = string(matches[1].first, matches[1].second);
		data->val = string(matches[2].first, matches[2].second);
		data->val_type = SCRIPT_VALUE::TYPE_NONE;

		pBlock->values.push_back(data);
		cp = (char*)matches[0].second;
	}
}

void SCRIPT_READER::ParseRegexScriptValue(const char* script, SCRIPT_VALUE::ValueType value_type, SCRIPT_BLOCK* pBlock)
{
	boost::regex* ex = 0;
	switch (value_type)
	{
	case SCRIPT_VALUE::TYPE_NUMBER:
		ex = const_cast<boost::regex*>(&REGEX_ASSIGN_NUMBER);
		break;
		
	case SCRIPT_VALUE::TYPE_VECTOR2:
		ex = const_cast<boost::regex*>(&REGEX_ASSIGN_VECTOR2);
		break;
	case SCRIPT_VALUE::TYPE_BOOLEAN:
		ex = const_cast<boost::regex*>(&REGEX_ASSIGN_BOOLEAN);
		break;
	case SCRIPT_VALUE::TYPE_STRING:
		ex = const_cast<boost::regex*>(&REGEX_ASSIGN_STRING);
		break;
	case SCRIPT_VALUE::TYPE_NONE:
		break;
	};

	char* cp = (char*)script;			 // Character point
	boost::cmatch matches;

	while (boost::regex_search(cp, matches, *ex))
	{
		SCRIPT_VALUE* data = new SCRIPT_VALUE();

		data->value_name = string(matches[1].first, matches[1].second);
		data->val = string(matches[2].first, matches[2].second);
		data->val_type = value_type;

		pBlock->values.push_back(data);
		cp = (char*)matches[0].second;
	}
}

bool SCRIPT_READER::SearchIncludedFile(const char* name)
{
	int i = 0;
	int size = (int)script_file_names.size();

	for (int i = 0; i < size; i++)
	{
		if (script_file_names[i].compare(name) == 0)
		{
			return true;
		}
	}

	return false;
}

void SCRIPT_READER::PreprocessIncludeSource(const char* src_script, string* dest_script)
{
	// Include, attach another source file
	boost::cmatch matches;

	string preprocessed_script = string(src_script);

	// Save the address for given data
	char* begin_script = (char*)preprocessed_script.data();

	PreprocessCommentSource(begin_script);		// Removing comment

	while (boost::regex_search(begin_script, matches, REGEX_INCLUDE))
	{
		string include_file_name = string(matches[1].first, matches[1].second);
		
		if (SearchIncludedFile(include_file_name.data()))
		{
			// if include_file_name is existed, delete source code
			char* begin_file_name = (char*)matches[0].first;
			char* end_file_name = (char*)matches[0].second;

			while (begin_file_name < end_file_name)
			{
				(*begin_file_name++) = ' ';
			}
			begin_file_name = (char*)matches[0].second;

			continue;
		}

		FILE* include_file;
		if (!(include_file = fopen(include_file_name.data(), "rt")))
		{
			char* begin_file_name = (char*)matches[0].first;
			char* end_file_name = (char*)matches[0].second;

			while (begin_file_name < end_file_name)
			{
				(*begin_file_name++) = ' ';
			}
			begin_script = (char*)matches[0].second;

			continue;
		}
		script_file_names.push_back(include_file_name);

		string include_script;
		include_script.clear();

		char char_include;
		while ((char_include = fgetc(include_file)) != EOF)
		{
			include_script += char_include;
		}

		string prev_script = string(preprocessed_script.data(), matches[0].first);
		string next_script = string(matches[0].second);

		preprocessed_script = prev_script+ '\n' + include_script + '\n' + next_script;
		PreprocessCommentSource(preprocessed_script.data());
	}

	(*dest_script) = preprocessed_script;
}

void SCRIPT_READER::PreprocessCommentSource(const char* script)
{
	char* cp = (char*)script;
	boost::cmatch matches;

	// Remove comment (// block name {), but regex is not stable
	boost::regex reg_comment_2("\\s*\\/\\/\\s*\\w*\\s*\\{\\s*");	// Note: \\s = spacing texture, \\/ = slash, \\w = alphabet including _, * = more than zero time
	
	while (boost::regex_search(cp, matches, reg_comment_2))
	{
		int block_count = 1;
		
		char* begin_comm = (char*)(matches[0].first);
		char* end_comm = (char*)(matches[0].second);

		while (block_count != 0)
		{
			if (*(end_comm) == '}')
			{
				block_count--;
			}
			if (*(end_comm) == '{')
			{
				block_count++;
			}
			end_comm++;
		}

		while (begin_comm < end_comm)
		{
			(*begin_comm++) = ' ';
		}
		cp = (char*)end_comm;
	}

	// Remove comment (//...), but regex is not stable
	boost::regex reg_comment_0("\\s*\\/\\/\\s*");
	cp = (char*)script;			// Character point
	while (boost::regex_search(cp, matches, reg_comment_0))
	{
		char* begin_comm = (char*)(matches[0].first);
		char* end_comm = (char*)(matches[0].second);

		while ((*end_comm) != '\0')
		{
			if ((*++end_comm) == '\n')
			{
				break;
			}
		}

		while (begin_comm < end_comm)
		{
			(*begin_comm++) = ' ';
		}
			
		cp = (char*)end_comm;
	}

	// Remove comment (/*...*/)
	boost::regex reg_comment_1("\\s*\\/\\*\\s*");
	cp = (char*)script;
	while (boost::regex_search(cp, matches, reg_comment_1))
	{
		char* begin_comm = (char*)(matches[0].first);
		char* end_comm = (char*)(matches[0].second);

		while ((*end_comm) != '\0')
		{
			if (*end_comm = '*' && *(end_comm + 1) == '/')
			{
				end_comm += 2;
				break;
			}
			end_comm++;
		}

		while (begin_comm < end_comm)
		{
			(*begin_comm++) = ' ';
		}
		cp = (char*)end_comm;
	}
}

void SCRIPT_READER::TrimSourceBlock(const char* script, vector<string>* trim_out)
{
	boost::regex reg("(\\w+)\\s*\\{\\s*");		// Note : \\w = alphabet including _, \\s = spacing character, * = more than zero time, + = more than one time 
	char* begin_syntax = (char*) script;
	char* end_syntax;
	boost::cmatch matches;

	int count = 0;

	while (boost::regex_search(begin_syntax, matches, reg))
	{
		begin_syntax = (char*)matches[0].first;
		end_syntax = (char*)matches[0].second;

		count = 1; // counting '{'
		while (count != 0)
		{
			if ((*end_syntax) == '{')
			{
				count++;
			}
			else if ((*end_syntax) == '}')
			{
				count--;
			}
			else if ((*end_syntax) == '\0')
			{
				break;
			}
			end_syntax++;
		}

		trim_out->push_back(string(begin_syntax, end_syntax));
		while (begin_syntax < end_syntax)
		{
			*(begin_syntax++) = ' ';
		}
	}
}

void SCRIPT_READER::TrimBlockNamePath(const char* path, vector<string>* trim_out)
{
	char* begin_path = (char*)path;
	char* end_path = (char*)path;

	while (*end_path != '\0')
	{
		if (*++end_path == '.')
		{
			trim_out->push_back(string(begin_path, end_path));
			begin_path = end_path + 1;
		}
	}

	trim_out->push_back(string(begin_path, end_path));
}

SCRIPT_BLOCK* SCRIPT_READER::SearchBlock(const vector<SCRIPT_BLOCK*>* blocks, const char* block_name) const
{
	if (!blocks)
	{
		return NULL;
	}

	int i = 0;
	int size = (int)blocks->size();

	for (int i = 0; i < size; i++)
	{
		if ((*blocks)[i]->block_name.compare(block_name) == 0)
		{
			return (*blocks)[i];
		}
	}

	return NULL;
}

SCRIPT_VALUE* SCRIPT_READER::SearchValue(const vector<SCRIPT_VALUE*>* values, const char* value_name) const
{
	if (!values)
	{
		return NULL;
	}

	int i = 0;
	int size = (int)values->size();

	for (int i = 0; i < size; i++)
	{
		if ((*values)[i]->value_name.compare(value_name) == 0)
		{
			return (*values)[i];
		}
	}

	return NULL;
}

const SCRIPT_BLOCK& SCRIPT_READER::FindBlock(const char* block_name)
{
	const SCRIPT_BLOCK* current_block = NULL;

	vector<string> trim;
	TrimBlockNamePath(block_name, &trim);
	
	if (current_block = SearchBlock(&script_blocks, trim[0].data()))
	{
		for (int i = 1; i < (int)trim.size(); i++)
		{
			current_block = SearchBlock(&(current_block->child_blocks), trim[i].data());
			if (!current_block)
			{
				return (*SCRIPT_READER::empty_block);
			}
		}
	}

	if (!current_block)
	{
		return (*SCRIPT_READER::empty_block);
	}
	
	return (*current_block);
}

SCRIPT_BLOCK* SCRIPT_READER::SearchBlock(const char* block_name)
{
	SCRIPT_BLOCK* current_block = NULL;

	vector<string> trim;
	TrimBlockNamePath(block_name, &trim);

	if (current_block = SearchBlock(&script_blocks, trim[0].data()))
	{
		for (int i = 1; i < (int)trim.size(); i++)
		{
			current_block = SearchBlock(&(current_block->child_blocks), trim[i].data());
			if (!current_block)
			{
				return 0;
			}
		}
	}

	if (!current_block)
	{
		return 0;
	}

	return current_block;
}

SCRIPT_BLOCK* SCRIPT_READER::AddNewBlock(const char* block_name)
{
	SCRIPT_BLOCK* new_block = new SCRIPT_BLOCK;
	new_block->block_name = string(block_name);

	script_blocks.push_back(new_block);
	return new_block;
}

SCRIPT_BLOCK* SCRIPT_READER::AddNewBlock(SCRIPT_BLOCK* block)
{
	script_blocks.push_back(block);
	return block;
}

SCRIPT_BLOCK* SCRIPT_READER::AddNewBlock(const char* block_name, const char* block_contents)
{
	SCRIPT_BLOCK* newBlock = AddNewBlock(block_name);
	newBlock->SetContents(block_contents);
	return newBlock;
}

void SCRIPT_READER::DeleteBlock(const char* block_name)
{
	vector<SCRIPT_BLOCK*>::iterator it = script_blocks.begin();
	for (it; it != script_blocks.end() ; it++)
	{
		if ((*it)->block_name.compare(block_name) == 0)
		{
			break;
		}
	}

	delete (*it);
	script_blocks.erase(it);
}

void SCRIPT_READER::WriteScript(const char* script_path)
{
	out_file_stream.open(script_path, ios_base::trunc | ios_base::out);
	assert(out_file_stream.is_open() == true);

	int size = (int)script_blocks.size();

	for (int i = 0; i < size; i++)
	{
		script_blocks[i]->WriteScript(0, out_file_stream);
	}

	out_file_stream.close();
}

SCRIPT_BLOCK* SCRIPT_READER::ParseSingleContents(const char* contents)
{
	ReleaseScriptReader();
	
	SCRIPT_BLOCK* new_block = new SCRIPT_BLOCK(); // creating new Data block
	script_blocks.push_back(new_block);

	ParseRegexScriptValue(contents, SCRIPT_VALUE::TYPE_NUMBER , new_block);
	ParseRegexScriptValue(contents, SCRIPT_VALUE::TYPE_VECTOR2, new_block);
	ParseRegexScriptValue(contents, SCRIPT_VALUE::TYPE_BOOLEAN, new_block);
	ParseRegexScriptValue(contents, SCRIPT_VALUE::TYPE_STRING, new_block);

	return new_block;		
}

void SCRIPT_BLOCK::TrimBlockNamePath(const char* path, vector<string>* trim_out) const
{
	char* begin_path = (char*)path;
	char* end_path = (char*)path;

	while (*end_path != '\0')
	{
		if (*++end_path == '.')
		{
			trim_out->push_back(string(begin_path, end_path));
			begin_path = end_path + 1;
		}
	}

	trim_out->push_back(string(begin_path, end_path));
}

const SCRIPT_BLOCK* SCRIPT_BLOCK::SearchBlock(const vector<SCRIPT_BLOCK*>* blocks, const char* block_name) const
{
	if (!blocks)
	{
		return NULL;
	}

	int i = 0;
	int size = (int)blocks->size();

	for (int i = 0; i < size; i++)
	{
		if ((*blocks)[i]->block_name.compare(block_name) == 0)
		{
			return (*blocks)[i];
		}
	}

	return NULL;
}

const SCRIPT_VALUE* SCRIPT_BLOCK::SearchValue(const vector<SCRIPT_VALUE*>* values, const char* value_name) const
{
	if (!values)
	{
		return NULL;
	}

	int i = 0;
	int size = (int)values->size();

	for (int i = 0; i < size; i++)
	{
		if ((*values)[i]->value_name.compare(value_name) == 0)
		{
			return (*values)[i];
		}
	}

	return NULL;
}

SCRIPT_BLOCK* SCRIPT_BLOCK::SearchBlock(const char* block_name)
{
	const SCRIPT_BLOCK* current_block = NULL;

	vector<string> trim;
	TrimBlockNamePath(block_name, &trim);

	if (current_block = SearchBlock(&child_blocks, trim[0].data()))
	{
		for (int i = 1; i < (int)trim.size(); i++)
		{
			current_block = SearchBlock(&(current_block->child_blocks), trim[i].data());
			if (!current_block)
			{
				return 0;
			}
		}
	}

	if (!current_block)
	{
		return 0;
	}

	return (SCRIPT_BLOCK*) current_block;
}

SCRIPT_VALUE* SCRIPT_BLOCK::SearchValue(const char* value_name)
{
	if (values.empty() == true)
	{
		return NULL;
	}

	int i = 0;
	int size = (int)values.size();

	for (int i = 0; i < size; i++)
	{
		if (values[i]->value_name.compare(value_name) == 0)
		{
			return values[i];
		}
	}

	return NULL;
}

SCRIPT_BLOCK* SCRIPT_BLOCK::AddNewBlock(const char* block_name)
{
	SCRIPT_BLOCK* new_block = new SCRIPT_BLOCK;
	new_block->block_name = string(block_name);

	child_blocks.push_back(new_block);
	return new_block;
}

SCRIPT_BLOCK* SCRIPT_BLOCK::AddNewBlock(SCRIPT_BLOCK* block)
{
	child_blocks.push_back(block);
	return block;
}

SCRIPT_BLOCK* SCRIPT_BLOCK::AddNewBlock(const char* block_name, const char* block_contents)
{
	SCRIPT_BLOCK* newBlock = AddNewBlock(block_name);
	newBlock->SetContents(block_contents);
	return newBlock;
}

SCRIPT_VALUE* SCRIPT_BLOCK::AddNewValue(const char* value_name, const char* value_string)
{
	SCRIPT_VALUE* new_value = new SCRIPT_VALUE;
	new_value->value_name = string(value_name);
	new_value->val = string(value_string);

	values.push_back(new_value);
	return new_value;
}

void SCRIPT_BLOCK::DeleteBlock(const char* block_name)
{
	vector<SCRIPT_BLOCK*>::iterator it = child_blocks.begin();
	for (it; it != child_blocks.end(); it++)
	{
		if ((*it)->block_name.compare(block_name) == 0)
		{
			break;
		}
	}

	delete (*it);
	child_blocks.erase(it);
}

void SCRIPT_BLOCK::DeleteValue(const char* value_name)
{
	vector<SCRIPT_VALUE*>::iterator it = values.begin();
	for (it; it != values.end(); it++)
	{
		if ((*it)->value_name.compare(block_name) == 0)
		{
			break;
		}
	}

	delete (*it);
	values.erase(it);
}
	
const SCRIPT_BLOCK& SCRIPT_BLOCK::FindBlock(const char* block_name) const
{
	// Set Data Block
	// How to use : BLOCK0.BLOCK1.BLOCK2
	const SCRIPT_BLOCK* current_block = NULL;

	vector<string> trim;
	TrimBlockNamePath(block_name, &trim);

	if (current_block = SearchBlock(&child_blocks, trim[0].data()))
	{
		for (int i = 1; i < (int)trim.size(); i++)
		{
			current_block = SearchBlock(&(current_block->child_blocks), trim[i].data());
			if (!current_block)
			{
				return (*SCRIPT_READER::empty_block);
			}
		}
	}

	if (!current_block)
	{
		return (*SCRIPT_READER::empty_block);
	}
		
	return (*current_block);
}

const SCRIPT_VALUE& SCRIPT_BLOCK::FindValue(const char* value_name) const
{
	if (values.empty() == true)
	{
		(*SCRIPT_READER::empty_value);
	}

	int size = (int)values.size();
	for (int i = 0; i < size; i++)
	{
		if (values[i]->value_name.compare(value_name) == 0)
		{
			return *(values[i]);
		}
	}

	return (*SCRIPT_READER::empty_value);
}

vector<SCRIPT_BLOCK*> SCRIPT_BLOCK::GetBlockList(const char* block_name)
{
	vector<SCRIPT_BLOCK*> blocks;
	const SCRIPT_BLOCK* current_block = 0;

	vector<string> trim;
	TrimBlockNamePath(block_name, &trim);

	int size = (int)child_blocks.size();
	for (int i = 0; i < size; i++)
	{
		if (child_blocks[i]->block_name.compare(block_name) == 0)
		{
			blocks.push_back(child_blocks[i]);
		}
	}

	return blocks;
}

int SCRIPT_BLOCK::GetInteger(const char* value_name, const int& default_value) const
{
	const SCRIPT_VALUE* data = SearchValue(&values, value_name);
	if (!data)
	{
		return default_value;
	}

	return atoi(data->val.data());
}

T SCRIPT_BLOCK::GetFloat(const char* value_name, const T& default_value) const
{
	const SCRIPT_VALUE* data = SearchValue(&values, value_name);
	if (!data)
	{
		return default_value;
	}

	return (T)atof(data->val.data());
}

VT SCRIPT_BLOCK::GetVector2(const char* value_name, const VT& default_value) const
{
	const SCRIPT_VALUE* data = SearchValue(&values, value_name);
	if (!data)
	{
		return default_value;
	}

	boost::cmatch matches;
	boost::regex_match(data->val.data(), matches, SCRIPT_READER::REGEX_ASSIGN_VECTOR2_VALUE);

	VT val;
		
	val.x = (T)atof(string(matches[1].first, matches[1].second).data());
	val.y = (T)atof(string(matches[2].first, matches[2].second).data());
	
	return val;
}



VI SCRIPT_BLOCK::GetInt2(const char* value_name, const VI& default_value) const
{
	const SCRIPT_VALUE* data = SearchValue(&values, value_name);
	if (!data)
	{
		return default_value;
	}

	boost::cmatch matches;
	boost::regex_match(data->val.data(), matches, SCRIPT_READER::REGEX_ASSIGN_VECTOR2_VALUE);

	VI val;

	val.x = atoi(string(matches[1].first, matches[1].second).data());
	val.y = atoi(string(matches[2].first, matches[2].second).data());
	
	return val;
}

bool SCRIPT_BLOCK::GetBoolean(const char* value_name, const bool default_value) const
{
	const SCRIPT_VALUE* data = SearchValue(&values, value_name);
	if (!data)
	{
		return default_value;
	}

	if (data->val.compare("true") == 0)
	{
		return true;
	}

	return false;
}

const char* SCRIPT_BLOCK::GetString(const char* value_name, const char* default_value) const
{
	const SCRIPT_VALUE* data = SearchValue(&values, value_name);
	if (!data)
	{
		return default_value;
	}

	return data->val.data();
}

void SCRIPT_BLOCK::WriteScript(const int level, ofstream& stream)
{
	for (int x = 0; x < level; x++)
	{
		stream << "\t";
		stream << block_name << endl;
	}

	for (int x = 0; x < level; x++)
	{
		stream << "\t";
		stream << "{" << endl;
	}

	int size = (int)values.size();

	for (int i = 0; i < size; i++)
	{
		for (int x = 0; x < level + 1; x++)
		{
			stream << "\t";
		}
		if (values[i]->val_type == SCRIPT_VALUE::TYPE_STRING)
		{
			stream << values[i]->value_name << " = " << "\"" << values[i]->val << "\"" << endl;
		}
		else
		{
			stream << values[i]->value_name << " = " << values[i]->val << endl;
		}
	}

	// TODO : Write sub blocks
	size = (int)child_blocks.size();
	for (int i = 0; i < size; i++)
	{
		child_blocks[i]->WriteScript(level + 1, stream);
	}
	for (int x = 0; x < level; x++)
	{
		stream << "\t";
		stream << "}" << endl;
	}
}

string SCRIPT_BLOCK::GetContents()
{
	string contents;
	int size = (int)values.size();
	for (int i = 0; i < size; i++)
	{
		string temp = values[i]->value_name;
		temp.append(" = ");
		if (values[i]->val_type == SCRIPT_VALUE::TYPE_STRING)
		{
			temp.append("\"");
		}
		temp.append(values[i]->val);
		if (values[i]->val_type == SCRIPT_VALUE::TYPE_STRING)
		{
			temp.append("\"");
		}
		if (i != (size - 1))
		{
			temp.append("\n");
		}
		contents.append(temp);
	}

	return contents;
}

void SCRIPT_BLOCK::SetContents(const string& contents)
{
	// Clear
	values.swap(vector<SCRIPT_VALUE*>());

	// Parse content
	SCRIPT_READER::ParseRegexScriptValue(contents.c_str(), SCRIPT_VALUE::TYPE_NUMBER, this);
	SCRIPT_READER::ParseRegexScriptValue(contents.c_str(), SCRIPT_VALUE::TYPE_VECTOR2, this);
	SCRIPT_READER::ParseRegexScriptValue(contents.c_str(), SCRIPT_VALUE::TYPE_BOOLEAN, this);
	SCRIPT_READER::ParseRegexScriptValue(contents.c_str(), SCRIPT_VALUE::TYPE_STRING, this);
}
