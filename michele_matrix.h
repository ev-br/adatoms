/*
 *  matrix.h
 *  
 *  thnx Michele Dolfi
 *
 */
#ifndef __MICHELE_MATRIX__
#define __MICHELE_MATRIX__
#include <vector>

class Matrix : public std::vector<double>
{
	size_t m;
	size_t n;
	
public:
	Matrix() : m(0), n(0) {}
	Matrix(size_t m_, size_t n_) : std::vector<double>(m_*n_), m(m_), n(n_) {}
	Matrix(size_t m_, size_t n_, double val) : std::vector<double>(m_*n_, val), m(m_), n(n_) {}
	
	void resize(size_t m_, size_t n_)
    {m = m_; n = n_; std::vector<double>::resize(m_*n_);}
	void reserve(size_t m_, size_t n_)
    {std::vector<double>::reserve(m_*n_);}
	void clear()
    {m = n = 0; std::vector<double>::clear();}
	
	size_t NRows() const {return m;}
	size_t NCols() const {return n;}
	
	double& operator()(size_t i, size_t j)
	{
		return operator[](i + j*m);
	}
	const double& operator()(size_t i, size_t j) const
	{
		return operator[](i + j*m);
	}
	void swap(Matrix &y)
	{
		std::vector<double>::swap(y);
		std::swap(n, y.n);
		std::swap(m, y.m);
	}
	void clearMemory()
	{
		Matrix empty;
		swap(empty);
	}
};

inline std::ostream& operator<< (std::ostream& os, const Matrix& M)
{
	for (size_t i=0; i<M.NRows(); ++i) {
		for (size_t j=0; j<M.NCols(); ++j) {
			os << M(i,j) << "\t";
		}
		os << std::endl;
	}
	return os;
}


#endif

