#pragma once

#include "COMMON_DEFINITION.h"
#include "VECTOR_ND.h"

#include <stdio.h>
#include <iostream>
#include <string>
#include <boost/regex.hpp>
using namespace std;

#define READ_SCRIPT_BLOCKS(script_block, itr_block)		vector<SCRIPT_BLOCK*>::iterator (itr_block);\
														for((itr_block) = (script_block).child_blocks.begin(); (itr_block) != (script_block).child_blocks.end();(itr_block)++)
class SCRIPT_BLOCK;

class SCRIPT_VALUE
{
public: // Enumerates
	enum ValueType
	{
		TYPE_NUMBER = 0,
		TYPE_VECTOR2,
		TYPE_BOOLEAN,
		TYPE_STRING,
		TYPE_NONE,
	};

public: // Essential Data
	string		value_name;
	string		val;
	ValueType	val_type;
};

class SCRIPT_READER
{
public: // Essential Data
	friend class SCRIPT_BLOCK;

	static const char*			REGEX_INCLUDE_STRING;
	static const char*			REGEX_BLOCKDATA_STRING;
	static const char*			REGEX_ASSIGN_NUMBER_STRING;
	static const char*			REGEX_ASSIGN_VECTOR2_STRING;
	static const char*			REGEX_ASSIGN_BOOLEAN_STRING;
	static const char*			REGEX_ASSIGN_VECTOR2_VALUE_STRING;
	static const char*			REGEX_ASSIGN_STRING_STRING;
	static const char*			REGEX_ASSIGN_FLAGS_STRING;
	static const char*			REGEX_FINDBLOCK_STRING;

	static const boost::regex	REGEX_FINDBLOCK;
	static const boost::regex	REGEX_INCLUDE;
	static const boost::regex	REGEX_BLOCKDATA;
	static const boost::regex	REGEX_VALUE_NAME;
	static const boost::regex	REGEX_ASSIGN_NUMBER;
	static const boost::regex	REGEX_ASSIGN_VECTOR2;
	static const boost::regex	REGEX_ASSIGN_BOOLEAN;
	static const boost::regex	REGEX_ASSIGN_STRING;
	static const boost::regex	REGEX_ASSIGN_VECTOR2_VALUE;
	static const boost::regex   REGEX_ASSIGN_FLAGS;

	static const SCRIPT_BLOCK*	empty_block;
	static const SCRIPT_VALUE*	empty_value;

	ofstream					out_file_stream;

	vector<SCRIPT_BLOCK*>		script_blocks;
	vector<string>				script_file_names;				// Included script file names

public: // Constructors and Destructor
	SCRIPT_READER(void)
	{}

	SCRIPT_READER(const char* script_filename)
	{
		Initialize(script_filename);
	}

	~SCRIPT_READER(void)
	{
		ReleaseScriptReader();
	}

public: // Initialization Function
	bool Initialize(const char* script);
	
public: // Member Functions
	void ReleaseScriptBlock(SCRIPT_BLOCK* block);
	void ReleaseScriptReader();
	
	void ParseRegexScriptBlock(const char* script, SCRIPT_BLOCK* pBlock = 0);
	static void ParseRegexScriptValue(const char* script, const boost::regex* ex, SCRIPT_BLOCK* pBlock = 0);
	static void ParseRegexScriptValue(const char* script, SCRIPT_VALUE::ValueType value_type, SCRIPT_BLOCK* pBlock = 0);
	
	bool SearchIncludedFile(const char* name);
	
	void PreprocessCommentSource(const char* script);
	void PreprocessIncludeSource(const char* src_script, string* dest_script);
	
	void TrimSourceBlock(const char* script, vector<string>* trim_out);
	void TrimBlockNamePath(const char* path, vector<string>* trim_out);
	
	SCRIPT_BLOCK* SearchBlock(const vector<SCRIPT_BLOCK*>* blocks, const char* block_name) const;
	SCRIPT_VALUE* SearchValue(const vector<SCRIPT_VALUE*>* values, const char* value_name) const;
	
	const SCRIPT_BLOCK& FindBlock(const char* block_name);
	vector<SCRIPT_BLOCK*>& GetBlocks();
	
	SCRIPT_BLOCK* SearchBlock(const char* block_name);
	SCRIPT_BLOCK* AddNewBlock(const char* block_name);
	SCRIPT_BLOCK* AddNewBlock(SCRIPT_BLOCK* block);
	SCRIPT_BLOCK* AddNewBlock(const char* block_name, const char* block_contents);
	void DeleteBlock(const char* block_name);
	
	void WriteScript(const char* script_path);
	
	SCRIPT_BLOCK* ParseSingleContents(const char* contents);
	
	friend class SCRIPT_BLOCK;
};

class SCRIPT_BLOCK
{
public: // Essential Data
	string						block_name;
	vector<SCRIPT_VALUE*>		values;
	vector<SCRIPT_BLOCK*>		child_blocks;

//public: // Constructor and Destructor
//	SCRIPT_BLOCK(void)
//	{}
//
//	~SCRIPT_BLOCK(void)
//	{}

public: // Member Functions
	void TrimBlockNamePath(const char* path, vector<string>* trim_out) const;
	
	const SCRIPT_BLOCK* SearchBlock(const vector<SCRIPT_BLOCK*>* blocks, const char* block_name) const;
	const SCRIPT_VALUE* SearchValue(const vector<SCRIPT_VALUE*>* values, const char* value_name) const;
	
	SCRIPT_BLOCK* SearchBlock(const char* block_name);
	SCRIPT_VALUE* SearchValue(const char* value_name);
	
	SCRIPT_BLOCK* AddNewBlock(const char* block_name);
	SCRIPT_BLOCK* AddNewBlock(SCRIPT_BLOCK* block);
	SCRIPT_BLOCK* AddNewBlock(const char* block_name, const char* block_contents);
	SCRIPT_VALUE* AddNewValue(const char* value_name, const char* value_string);
	
	void DeleteBlock(const char* block_name);
	void DeleteValue(const char* value_name);
		
	const SCRIPT_BLOCK& FindBlock(const char* block_name) const;
	const SCRIPT_VALUE& FindValue(const char* value_name) const;
	
	vector<SCRIPT_BLOCK*> GetBlockList(const char* block_name);
	
	int				GetInteger(const char* value_name, const int& default_value = 0							) const;
	T				GetFloat(const char* value_name, const T& default_value = (T)0							) const;
	VT				GetVector2(const char* value_name, const VT& default_value = VT()						) const;
	VI				GetInt2(const char* value_name, const VI& default_value = VI()							) const;
	bool			GetBoolean(const char* value_name, const bool default_value = false						) const;
	const char*		GetString(const char* value_name, const char* default_value = "Null"					) const;
	
	int GetValue(const char* value_name, const int& default_value = 0) const
	{
		return GetInteger(value_name, default_value);
	}

	T GetValue(const char* value_name, const T& default_value = (T)0) const
	{
		return GetFloat(value_name, default_value);
	}

	VT GetValue(const char* value_name, const VT& default_value = VT()) const
	{
		return GetVector2(value_name, default_value);
	}

	VI GetValue(const char* value_name, const VI& default_value = VI()) const
	{
		return GetInt2(value_name, default_value);
	}

	bool GetValue(const char* value_name, const bool default_value = false) const
	{
		return GetBoolean(value_name, default_value);
	}

	const char* GetValue(const char* value_name, const char* default_value = "Null") const
	{
		return GetString(value_name, default_value);
	}

	void WriteScript(const int level, ofstream& stream);
	
	string GetContents();
	void SetContents(const string& contents);
};

static void FileText(const char* file, string& text)
{
	FILE* file_script = 0;
	file_script = fopen(file, "rt");

	assert(file != NULL);
	text.clear();

	char char_script;
	while ((char_script = fgetc(file_script)) != EOF)
	{
		text += char_script;
	}

	fclose(file_script);
}


	


		