#pragma once

#include "COMMON_DEFINITION.h"
#include "LEVELSET_OBJECT.h"

class PLANE : public LEVELSET_OBJECT
{
public: // Essential Data
	VT position;
	VT unit_normal;

public: // Constructors and Destructor
	PLANE(void)
	{}

	PLANE(const VT& position_input, const VT& normal_input)
	{
		Initialize(position_input, normal_input);
	}
	
	~PLANE(void)
	{}

public: // Initialization Function
	void Initialize(const VT& position_input, const VT& normal_input)
	{
		position = position_input;
		unit_normal = normal_input;
		unit_normal.MakeThisUnit();
	}

public: // Member Functions
	const T SignedDistance(const VT& sampling_location) const
	{
		return DotProduct(sampling_location - position, unit_normal);
	}

	void Normal(const VT& position, VT& normal) const
	{
		normal = unit_normal;
	}

	const VT Normal(const VT& position) const
	{
		return unit_normal;
	}

	void UnitNormal(const VT& sampling_location, VT& normal_output) const
	{
		normal_output = unit_normal;
	}

	const VT UnitNormal(const VT& position) const
	{
		return unit_normal;
	}
};






