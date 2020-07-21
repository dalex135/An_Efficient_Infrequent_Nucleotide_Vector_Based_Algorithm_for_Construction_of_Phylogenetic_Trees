#pragma once

#include <iostream>
#include "DNA_UTILS.h"

class btree
{

private:
	struct node
	{
		DNA_UTILS::differ key_value;
		node *left;
		node *right;
	};
	void destroy_tree(node *leaf);
	void insert(DNA_UTILS::differ key, node *leaf);
	node *search(int key, node *leaf);
	void outTree(node *);
	node *root;
	vector<DNA_UTILS::differ> sortedTree;
public:
	btree();
	~btree();

	void insert(DNA_UTILS::differ key);
	node *search(int key);
	void destroy_tree();
	vector<DNA_UTILS::differ> printTree();
};

