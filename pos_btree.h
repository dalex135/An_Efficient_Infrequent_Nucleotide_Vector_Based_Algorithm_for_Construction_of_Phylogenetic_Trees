#pragma once
#include <iostream>
#include "DNA_UTILS.h"

class pos_btree
{
	struct node
	{
		DNA_UTILS::positiondata key_value;
		node *left;
		node *right;
	};
	void destroy_tree(node *leaf);
	void insert(DNA_UTILS::positiondata key, node *leaf);
	node *search(DNA_UTILS::positiondata key, node *leaf);
	void outTree(node *);
	node *root;
	vector<DNA_UTILS::positiondata> sortedTree;
public:
	pos_btree();
	~pos_btree();

	void insert(DNA_UTILS::positiondata key);
	node *search(DNA_UTILS::positiondata key);
	void destroy_tree();
	vector<DNA_UTILS::positiondata> printTree();
};

