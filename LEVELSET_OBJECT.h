#pragma once

#include "COMMON_DEFINITION.h"

// Need to know why they define this class
// Abstract class for BOX, SPHERE, CYLINDER, PLANE
class LEVELSET_OBJECT
{
public: // Constructors and Destructor
	LEVELSET_OBJECT(void)
	{}

	~LEVELSET_OBJECT(void)
	{}

public: // Operator Overloading
	const T operator()(const VT& position) const
	{
		return SignedDistance(position);
	}
	
	const T operator()(const int& i, const int& j) const
	{
		return SignedDistance(i, j);
	}

public: // Virtual Functions for abstract class of BOX, SPHERE, CYLINDER, PLANE
	virtual const T SignedDistance(const VT& position) const
	{
		cout << "Virtual Function LEVELSET_OBJECT::SignedDistace(const VT& position) was called" << endl;
		exit(-1);
		return T();
	}

	virtual const T SignedDistance(const int& i, const int& j) const
	{
		cout << "Virtual Function LEVELSET_OBJECT::SignedDistace(const int& i, const int& j) was called" << endl;
		exit(-1);
		return T();
	}
	virtual void Normal(const VT& position, VT& normal) const
	{
		cout << "Virtual Function LEVELSET_OBJECT::Normal(const VT& position, VT& normal) was called" << endl;
		exit(-1);
	}

	virtual const VT Normal(const VT& position) const
	{
		cout << "Virtual Function LEVELSET_OBJECT::Normal(const VT& position) was called" << endl;
		exit(-1);
		return VT();
	}
		
	virtual void UnitNormal(const VT& position, VT& normal) const
	{
		cout << "Virtual Function LEVELSET_OBJECT::UnitNormal(const VT& position, VT& normal) was called" << endl;
		exit(-1);
	}

	virtual const VT UnitNormal(const VT& position) const
	{
		cout << "Virtual Function LEVELSET_OBJECT::UnitNormal(const VT& position) was called" << endl;
		exit(-1);
		return VT();
	}

	virtual void Update() const // Recompute Normal array here
	{
		cout << "Virtual Function LEVELSET_OBJECT::Update() was called" << endl;
		exit(-1);
	}

	void Surface(const VT& position, VT& surface_position) const
	{
		surface_position = position;
		VT unitnormal(UnitNormal(position));
		T phi = SignedDistance(position);
		surface_position -= unitnormal*phi;
	}
};



	

