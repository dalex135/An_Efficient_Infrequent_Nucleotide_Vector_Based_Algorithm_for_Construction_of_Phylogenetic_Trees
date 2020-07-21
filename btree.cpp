#include "btree.h"


btree::btree()
{
	root = NULL;
}


btree::~btree()
{
	sortedTree.clear();
	destroy_tree();
}

void btree::destroy_tree(node *leaf)
{
	if (leaf != NULL)
	{
		destroy_tree(leaf->left);
		destroy_tree(leaf->right);
		delete leaf;
	}
}

void btree::insert(DNA_UTILS::differ key, node *leaf)
{
	if (key.gap == leaf->key_value.gap)
	{
		leaf->key_value.freaqency++;
		leaf->key_value.direct = leaf->key_value.direct + key.direct;
		leaf->key_value.neighbor = leaf->key_value.neighbor + key.neighbor;
	}
	else if (key.gap < leaf->key_value.gap)
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
	else if (key.gap > leaf->key_value.gap)
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

btree::node *btree::search(int key, node *leaf)
{
	if (leaf != NULL)
	{
		if (key == leaf->key_value.gap)
			return leaf;
		if (key < leaf->key_value.gap)
			return search(key, leaf->left);
		else
			return search(key, leaf->right);
	}
	else return NULL;
}

void btree::insert(DNA_UTILS::differ key)
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

void btree::outTree(node *track)
{
	if (track == NULL)
		return;
	outTree(track->left);
	sortedTree.push_back(track->key_value);
	outTree(track->right);
}

vector<DNA_UTILS::differ> btree::printTree()
{
	sortedTree.clear();
	outTree(root);
	return sortedTree;
}

void btree::destroy_tree()
{
	destroy_tree(root);
}