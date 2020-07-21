#include "pos_btree.h"


pos_btree::pos_btree()
{
	root = NULL;
}


pos_btree::~pos_btree()
{
	sortedTree.clear();
	destroy_tree();
}

void pos_btree::destroy_tree(node *leaf)
{
	if (leaf != NULL)
	{
		destroy_tree(leaf->left);
		destroy_tree(leaf->right);
		delete leaf;
	}
}

void pos_btree::insert(DNA_UTILS::positiondata key, node *leaf)
{
	if (key.a_pos == leaf->key_value.a_pos && key.b_pos == leaf->key_value.b_pos)
		return;
	else if (key.a_pos < leaf->key_value.a_pos)
	{
		if (leaf->left != NULL)
			insert(key, leaf->left);
		else
		{
			leaf->left = new node;
			leaf->left->key_value = key;
			leaf->left->left = NULL;    //Sets the left child of the child node to null
			leaf->left->right = NULL;   //Sets the right child of the child node to null
		}
	}
	else if (key.a_pos > leaf->key_value.a_pos)
	{
		if (leaf->right != NULL)
			insert(key, leaf->right);
		else
		{
			leaf->right = new node;
			leaf->right->key_value = key;
			leaf->right->left = NULL;  //Sets the left child of the child node to null
			leaf->right->right = NULL; //Sets the right child of the child node to null
		}
	}
}

pos_btree::node *pos_btree::search(DNA_UTILS::positiondata key, node *leaf)
{
	if (leaf != NULL)
	{
		if (key.a_pos == leaf->key_value.a_pos && key.b_pos == leaf->key_value.b_pos)
			return leaf;
		if (key.a_pos < leaf->key_value.a_pos)
			return search(key, leaf->left);
		else
			return search(key, leaf->right);
	}
	else return NULL;
}

void pos_btree::insert(DNA_UTILS::positiondata key)
{
	if (root != NULL)
		insert(key, root);
	else
	{
		root = new node;
		root->key_value = key;
		root->left = NULL;
		root->right = NULL;
	}
}

void pos_btree::outTree(node *track)
{
	if (track == NULL)
		return;
	outTree(track->left);
	sortedTree.push_back(track->key_value);
	outTree(track->right);
}

vector<DNA_UTILS::positiondata> pos_btree::printTree()
{
	sortedTree.clear();
	outTree(root);
	return sortedTree;
}

void pos_btree::destroy_tree()
{
	destroy_tree(root);
}
