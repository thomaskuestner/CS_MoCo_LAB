#include "SubsampleMain.h"


LinkedList::LinkedList(void)
{
}
LinkedList::~LinkedList(void)
{
}

void LinkedList::insertElement(Point newPoint, LinkedList *&anchor, int &nElements)
{
	// list is empty
	if(anchor == 0)
	{
		LinkedList* element = new LinkedList();
		element->data = newPoint;
		element->next = 0;
		anchor = element;
	}
	// list is not empty
	else
	{
		LinkedList* element = anchor;
		LinkedList* newElement = new LinkedList();
		while(element->next != 0)
			element = element->next;
		newElement->data = newPoint;
		newElement->next = 0;
		element->next = newElement;
	}
	nElements++;
}

void LinkedList::deleteElement(int position, LinkedList *&anchor, int &nElements)
{
	if(position < nElements)
	{
		// delete first element in list
		if(position == 0)
		{
			LinkedList* del = anchor;
			anchor = anchor->next;
			delete del;
		}
		// iterate to desired position and delete element
		else
		{
			LinkedList* thisElement = anchor;
			LinkedList* elementBefore = anchor;
			for(int i=0; i<position; i++)
			{
				elementBefore = thisElement;
				thisElement = thisElement->next;
			}
			LinkedList* elementBehind = thisElement->next;
			elementBefore->next = elementBehind;
			delete thisElement;
		}
		nElements--;
	}
	else
	{
		cout << "error" << endl << endl;
		system("PAUSE");
		exit(1);
	}
}

Point LinkedList::showElement(int position, LinkedList *anchor, int nElements)
{
	if(position < nElements)
	{
		Point showThis;
		LinkedList* thisElement = anchor;
		LinkedList* element = new LinkedList();

		for(int i=0; i<position; i++)
		{
			thisElement = thisElement->next;
		}
		if(thisElement != 0)
			showThis = thisElement->data;

		return showThis;
	}
	else
	{
		cout << "error" << endl << endl;
		system("PAUSE");
		exit(1);
	}
}
