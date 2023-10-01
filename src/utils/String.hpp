﻿#pragma once
#define __STDC_WANT_LIB_EXT1__ 1
#include <io.h>
#include <string>
#include <string.h>

namespace str_util
{
	using std::string;

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#  define DELIMITER "\\"
#else
#  define DELIMITER "/"
#endif

#ifndef CONCACT(x,y)
#  define CONCACT(x,y) x##y
#endif

	template<typename T = string>
	string concatFilePath(T firstArg)
	{
		return firstArg;
	}

	// 目前只支持参数全为string类型的情况，待处理参数包中携带const char*(const char[])的这种情况
	template<typename T = string, typename...Types>
	string concatFilePath(T firstArg, Types ... args)
	{
		firstArg = firstArg + (string)DELIMITER;
		return firstArg + concatFilePath(args...);
	}

	inline char* getFileNameWithExt(const char delimiter, char* filePath)
	{
		char* fileName = filePath;
		const size_t len = strlen(filePath);
		for (int i = 0; i < len; ++i)
		{
			if (filePath[i] == delimiter)
				fileName = &filePath[i + 1];
		}
		return fileName;
	}

	inline char* getFileNameWithExt(const char* delimiter, const char* filePath)
	{
		char* filePath_c = new char[strlen(filePath) + 1];
		strcpy(filePath_c, filePath);
		char* fileName = nullptr;

#ifdef __STDC_LIB_EXT1__
		char* token = strtok(filePath_c, &delimiter); // 'filePath_c' is changed after calles strtok
		while (token != nullptr)
		{
			if (fileName != nullptr) { delete[] fileName; fileName = nullptr; }
			fileName = new char[strlen(token) + 1];
			strcpy(fileName, token);
			token = strtok(NULL, &delimiter); // strtok第一个参数传入nullptr，代表使用之前strtok内部保存的SAVE_PTR定位到下一个待处理的字符的位置
		}
#else // Thread safe
		char* ptr = nullptr;
		char* token = strtok_s(filePath_c, delimiter, &ptr); //相较于strtok()函数，strtok_s函数需要用户传入一个指针，用于函数内部判断从哪里开始处理字符串
		while (token != nullptr)
		{
			if (fileName != nullptr) { delete[] fileName; fileName = nullptr; }
			fileName = new char[strlen(token) + 1];
			strcpy(fileName, token);
			token = strtok_s(NULL, delimiter, &ptr);
		}
#endif
		delete[] filePath_c;
		filePath_c = nullptr;
		return fileName;
	}

	inline string getFileNameWithExt(const char* delimiter, const string& filePath)
	{
		string fileName = getFileNameWithExt(delimiter, filePath.c_str());
		return fileName;
	}

	// such as ".obj"
	inline char* getFileExtension(char* filePath)
	{
		char* extension = filePath;
		for (int i = strlen(filePath); i >= 0; --i)
		{
			if (filePath[i] == '.')
			{
				extension = &filePath[i + 1];
				return extension;
			}
		}
		extension[0] = 0;
		return extension;
	}

	inline char* getFileExtension(const char* filePath)
	{
		char* temp = new char[strlen(filePath) + 1];
		strcpy(temp, filePath);
		char* extension = getFileExtension(temp);
		delete[] temp;
		temp = nullptr;
		return extension;
	}

	inline string getFileExtension(const string& filePath)
	{
		const size_t idx = filePath.rfind('.');
		if (idx != string::npos)
			return filePath.substr(idx);
		return "";
	}

	inline char* getFileName(const char* delimiter, char* filePath)
	{
		char* fileName = getFileNameWithExt(delimiter, filePath);
		size_t len = strlen(fileName);
		for (int i = len - 1; i >= 0; --i)
		{
			if (fileName[i] == '.')
			{
				fileName[i] = 0;
				return fileName;
			}
		}
		return fileName;
	}

	inline char* getFileName(const char* delimiter, const char* filePath)
	{
		char* fileName = getFileNameWithExt(delimiter, filePath);
		size_t len = strlen(fileName);
		for (int i = len - 1; i >= 0; --i)
		{
			if (fileName[i] == '.')
			{
				fileName[i] = 0;
				return fileName;
			}
		}
		return fileName;
	}

	inline string getFileName(const char* delimiter, const string& filePath)
	{
		string fileName = getFileName(delimiter, filePath.c_str());
		return fileName;
	}

	inline char* getDirName(const char* delimiter, const char* filePath)
	{
		const size_t len = strlen(filePath);
		char* dirName = new char[len + 1];
		strcpy(dirName, filePath);
		for (int i = len - 1; i >= 0; --i)
		{
			if (filePath[i] == delimiter[0])
			{
				dirName[i] = 0;
				return dirName;
			}
		}
		dirName[0] = 0;
		return dirName;
	}

	inline char* getDirName(const char* delimiter, char* filePath)
	{
		const size_t len = strlen(filePath);
		char* dirName = new char[len + 1];
		strcpy(dirName, filePath);
		for (int i = len - 1; i >= 0; --i)
		{
			if (filePath[i] == delimiter[0])
			{
				dirName[i] = 0;
				return dirName;
			}
		}
		dirName[0] = 0;
		return dirName;
	}

	inline void checkDir(const string& filename)
	{
		size_t dir_idx = filename.find_last_of("\\");
		if (dir_idx != string::npos)
		{
			//cout << "output dir = " << filename.substr(0, dir_idx) << endl;
			string dir = filename.substr(0, dir_idx);
			if (_access(dir.c_str(), 0) == -1)
			{
				string command = "mkdir " + dir;
				system(command.c_str());
			}
		}
	}
}// namespace str_util