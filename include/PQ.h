
#ifndef PriorityQueue_H
#define PriorityQueue_H


typedef unsigned int _index;
typedef unsigned int integer;
typedef float  sreal;
typedef double dreal;

// specific data structures
template <class Element>
class PriorityQueue
{
protected:
	typedef std::map<Element*, _index> access_t;
	Element** S;
	access_t I;
	_index size;
	_index maxSize;
public:
	PriorityQueue() {}
	PriorityQueue(_index n);

	void swap(_index i, _index j);
	bool isEmpty() const;

	_index left(_index i) const;
	_index right(_index i) const;
	_index parent(_index i) const;


	void minHeapify(_index i);
	void decreaseKey(_index ind, sreal value);
	void increaseKey(_index ind, sreal value);
	void changeKey(Element* elem, sreal newvalue);
	void deleteKey(_index ind);
	void deleteElement(Element* e);
	void push(Element* elem);
	void pop();
	//Element* operator[]{}
	


	Element* top() const;
	//void print() const;

	~PriorityQueue();

};

template <class Element> void PriorityQueue<Element>::swap(_index i, _index j)
{
	Element* e_temp = S[i];
	_index i_temp = I[e_temp];
	I[S[i]] = I[S[j]];
	S[i] = S[j];
	S[j] = e_temp;
	I[S[i]] = i_temp;
}
// implementation
template <class Element> PriorityQueue<Element>::PriorityQueue(_index n) : size(0), maxSize(n)
{
	S = new Element*[n];
}

template <class Element> bool PriorityQueue<Element>::isEmpty() const
{
	return this->size==0;
}

template <class Element> _index PriorityQueue<Element>::left(_index i) const
{
	_index leftChild = 2 * i + 1;
	if (leftChild >= this->size)  return i;
	return leftChild;
}

template <class Element> _index PriorityQueue<Element>::right(_index i) const
{
	_index rightChild = 2 * i + 2;
	if (rightChild >= this->size) return i;
	return rightChild;
}

template <class Element> _index PriorityQueue<Element>::parent(_index i) const
{
	return (i-1) / 2;
}

template <class Element> void PriorityQueue<Element>::minHeapify(_index i)
{
	_index smallest, l, r;
	while (true)
	{
		l = left(i);
		r = right(i);
		if (*S[l]<*S[i]) { smallest = l; }
		else { smallest = i; }
		if (*S[r]<*S[smallest]) smallest = r;
		if (smallest == i) break;
		// exchange elements in container
		swap(i,smallest);
		i = smallest;
	}
	return;
}


template <class Element> void PriorityQueue<Element>::decreaseKey(_index ind, sreal value)
{
	if (value >= S[ind]->get_price()) return;
	S[ind]->set_price(value);
	while (ind>0 && *S[ind]<*S[parent(ind)])
	{
		swap(ind, parent(ind));
		ind = parent(ind);
	}
	return;
}

template <class Element> void PriorityQueue<Element>::increaseKey(_index ind, sreal value)
{
	if (value <= S[ind]->get_price()) return;
	S[ind]->set_price(value);
	minHeapify(ind);
	return;
}

template <class Element>void PriorityQueue<Element>::changeKey(Element* elem, sreal newvalue)
{
	if (newvalue > S[I[elem]]->get_price()) increaseKey(I[elem], newvalue);
	if (newvalue < S[I[elem]]->get_price()) decreaseKey(I[elem], newvalue);
	return;
}

template <class Element> void PriorityQueue<Element>::deleteKey(_index ind)
{
	if (ind >= this->size) return;
	I.erase(S[ind]);
	S[ind] = S[this->size - 1];
	I[S[ind]] = ind;
	this->size = this->size - 1;
	minHeapify(ind);
	return;
}

template <class Element> void PriorityQueue<Element>::deleteElement(Element* e)
{
	if (!I.count(e)) return;
	deleteKey(I[e]);
	return;
}
template <class Element> void PriorityQueue<Element>::push(Element* elem)
{
	this->size = this->size + 1;
	dreal temp = elem->get_price();
	elem->set_price(std::numeric_limits<dreal>::infinity());
	S[this->size - 1] = elem;
	I[elem] = this->size - 1;
	decreaseKey(this->size - 1, (sreal)temp);
	return;
}

template <class Element> void PriorityQueue<Element>::pop()
{
	if (this->size<1) { std::cout << "Priority Queue is empty!" << std::endl; return; }
	I.erase(S[0]);
	S[0] = S[this->size - 1];
	I[S[0]] = 0;
	this->size = this->size - 1;
	minHeapify(0);
	return;
}

template <class Element> Element* PriorityQueue<Element>::top() const
{
	if (this->size > 0) return S[0];
	std::cout << "\nPriority Queue is empty. NULL is returned." << std::endl;
	return NULL;
}

/*template <class Element> void PriorityQueue<Element>::print() const
{
for (_index i = 0; i<this->size; i++) S[i]->print();
std::cout << endl;
}*/




template <class Element> PriorityQueue<Element>::~PriorityQueue()
{
	for (_index i = 0; i<this->size; i++)
	{
		delete S[i];
	}
	delete S;
}

#endif