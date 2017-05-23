class Particle
{
	public:
	int type;
	Point3D getPosition();				//orizetai Point3D dioti epistrefei ena antikeimeno Point3D
	void setPosition();
	std::vector getClosestNeighbors();	
	
	
	private:
	Point3D p;
	int set;
	std::vector<int> closest_neighbors(9,0);
};