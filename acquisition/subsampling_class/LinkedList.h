#ifndef LINKEDLIST_H
#define LINKEDLIST_H

class LinkedList
{
public:
	LinkedList(void);
	virtual ~LinkedList(void);

	Point data;
	LinkedList* next;

	static void insertElement(Point newPoint, LinkedList *&anchor, int &nElements);
	static void deleteElement(int position, LinkedList *&anchor, int &nElements);
	static Point showElement(int position, LinkedList *anchor, int nElements);
};


#endif