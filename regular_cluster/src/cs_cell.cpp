#include "cs_cell.hpp"

//---------------------------------------------------------------------------------------------------------------------------------

Cell::Cell( glm::dvec3 orig, double radius )
{
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> distributionID( 1,maxID );

    origin = orig;
    R = radius;
    maxElongation = 0.7*R;
    name = distributionID(generator);
    age = 0.0;
    conf_num = 0;
    isLinked = false;
    linkedID = 0;
    isElongating = false;
    isVirtual = false;
    currentElongation = 0.0;
    pastElongation = 0.0;
    growthRate = growthRateConstant;
    conf_degree = 0.0;
    conf_degree_upd = 0.0;

    x = origin.x;
    y = origin.y;
    z = origin.z;

}

//---------------------------------------------------------------------------------------------------------------------------------

Cell::~Cell()
{
}

//---------------------------------------------------------------------------------------------------------------------------------

Cell Cell::addcell()
{
    Cell c( origin_nov,R );
    c.isLinked = true;
    c.linkedID = name;
    c.isElongating = true;
    c.growthdir = growthdir;
    c.growthRate = growthRate_nov;
    c.age = age;
    c.isVirtual = true;
    c.conf_num = conf_num;
    c.conf_degree = conf_degree;
    c.conf_degree_upd = conf_degree_upd;
    //c.confinements = confinements;
    isLinked = true;
    linkedID = c.name;

    c.Ath = Ath;
    c.Ax = Ax;
    c.Ay = Ay;
    c.Az = Az;
    c.Fth = Fth;
    c.Fx = Fx;
    c.Fy = Fy;
    c.Fz = Fz;
    c.p = p;
    c.th = th;
    c.Vth = Vth;
    c.Vx = Vx;
    c.Vy = Vy;
    c.Vz = Vz;
    c.normOfSumF = normOfSumF;
    c.sumOfNormF = sumOfNormF;
    c.isActive = isActive;
    c.m = m;
    c.rhoInd = rhoInd;
    c.I = I;

    return c;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Cell::growth( std::vector<Cell> &cells, double time )
{
    age = age + time;

    conf_degree_upd = updConfDegree( cells );

    if ( isElongating ) {

        if ( !isLinked ) return;

        if ( currentElongation/maxElongation >= 1.0 ) {
            /* release daughter cell */
            pastElongation = currentElongation;
            isLinked = false;
            isVirtual = false;
            linkedID = 0;
            age = 0.0;
            //conf_num = 0;
            //conf_degree = 0.0;
            isElongating = false;
            currentElongation = 0.0;
            growthRate = growthRateConstant;
            confinements.clear();
            return;
        }

        origin = getNewOrigin( origin, growthdir, growthRate );
        x = origin.x;
        y = origin.y;
        z = origin.z;
        for (size_t i = 0; i < cells.size(); i++) {
            if ( name == cells[i].linkedID ) {
                if ( isVirtual ){
                    currentElongation = cells[i].currentElongation;
                    break;
                }
                glm::dvec3 v1 = origin;
                glm::dvec3 v2 = cells[i].origin;
                glm::dvec3 v12 = v1 - v2;
                currentElongation = glm::length(v12);
                break;
            }
        }
    }
    else {
        if ( age >= initElongationTime ) {

            getGrowthDirections( cells );

            if ( conf_degree > growthThreshold ){
                conf_degree = 1.0;
            }

            switch ( conf_num )
            {
            case 0: case 1: case 2:
                growthRate = growthRateConstant * ( 1.0 - conf_degree );
                growthRate_nov = -growthRateConstant * ( 1.0 - conf_degree );
                break;        
            default:
                growthRate = 0.0;
                growthRate_nov = -2.0 * growthRateConstant * ( 1.0 - conf_degree );
                break;
            }

            divisiondir = growthdir;
            glm::dvec3 currorigin = origin;
            origin = getNewOrigin( currorigin, growthdir, growthRate );
            origin_nov = getNewOrigin( currorigin, growthdir, growthRate_nov );
            isElongating = true;
            x = origin.x;
            y = origin.y;
            z = origin.z;
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------------------

glm::dvec3 Cell::getNormale( glm::dvec3 p0, glm::dvec3 p1, glm::dvec3 p2 )
{
    if ( glm::length(p0) == 0.0 ) {
        p0.x = p0.x + 0.0001;
        p0.y = p0.y + 0.0001;
        p0.z = p0.z + 0.0001;
    }

    glm::dvec3 v1 = p1 - p0;
    glm::dvec3 v2 = p2 - p0;
    glm::dvec3 normale = glm::cross( v1, v2 );

    if ( glm::length(normale) == 0.0 )
        return normale;
    else
        return glm::normalize(normale);       
}

//---------------------------------------------------------------------------------------------------------------------------------

glm::dvec3 Cell::getCross( glm::dvec3 p1, glm::dvec3 p2 )
{
    glm::dvec3 v = glm::cross( p1, p2 );

    if ( glm::length(v) == 0.0 )
        return v;
    else
        return glm::normalize(v);       
}

//---------------------------------------------------------------------------------------------------------------------------------

glm::dvec3 Cell::getOrthVect( glm::dvec3 p )
{
    if ( p.x == 0.0 )
        p.x = 0.00001;
    if ( p.y == 0.0 )
        p.y = 0.00001;
    return glm::dvec3(1.0/p.x, 1.0/p.y, 0.0);
}

//---------------------------------------------------------------------------------------------------------------------------------

void Cell::getConfinements( std::vector<Cell> &cells )
{
    std::vector<glm::dvec3> t_confinements;
    bool isModified = false;

    for ( size_t j = 0; j < cells.size(); j++ ){
        if ( name != cells[j].name ){
        /*  more precise condition could be:
            if virtual cell is mature enough, then
            take it into account for confinement */
            // if ( cells[j].isVirtual && (cells[j].currentElongation/cells[j].maxElongation < 0.5) ){ /* a cell is elongated enough to make impact on confinement */
            //     continue;
            // }
            glm::dvec3 v1 = origin;
            glm::dvec3 v2 = cells[j].origin;
            glm::dvec3 v12 = v1 - v2;
            double d = glm::length(v12); //std::cout<<"d = "<<d/R<<std::endl;
            if (  d < R*confLimit ) { /* can be called 'response to confinement', probably can be responsible for an aggregate shape */
                glm::dvec3 v;
                v.x = v12.x/d;
                v.y = v12.y/d;
                v.z = v12.z/d;
                t_confinements.push_back( v );
                isModified = true;
            }
        }
    }

    if (isModified){
        size_t tc_size = t_confinements.size();
        if ( tc_size > 2 ){
            for (size_t i = 0; i < tc_size; i++){
                for (size_t j = 0; j < tc_size; j++){
                    if ( i != j ){
                        double res = std::abs( glm::dot( t_confinements[i],t_confinements[j] ) );
                        if ( res < 0.05 ){
                            confinements.push_back(t_confinements[i]);
                            break;
                        }
                    }
                }
            }
        }
        else {
            confinements = t_confinements;
        }
    }

    conf_num = confinements.size();
    //if (cells.size() == 4) std::cout<<"conf_num = "<<conf_num<<std::endl;
    conf_degree = 1.0 - 1.0 / ( (double)(conf_num + 1) );

    t_confinements.clear();
    t_confinements.shrink_to_fit();
}

//---------------------------------------------------------------------------------------------------------------------------------

double Cell::updConfDegree( std::vector<Cell> &cells )
{
    std::vector<glm::dvec3> t_confinements;
    bool isModified = false;

    for ( size_t j = 0; j < cells.size(); j++ ){
        if ( name != cells[j].name ){
        /*  more precise condition could be:
            if virtual cell is mature enough, then
            take it into account for confinement */
            if ( cells[j].isVirtual && (cells[j].currentElongation/cells[j].maxElongation < 0.5) ){ /* a cell is elongated enough to make impact on confinement */
                continue;
            }
            glm::dvec3 v1 = origin;
            glm::dvec3 v2 = cells[j].origin;
            glm::dvec3 v12 = v1 - v2;
            double d = glm::length(v12); //std::cout<<"d = "<<d/R<<std::endl;
            if (  d < R*confLimit ) { /* can be called 'response to confinement', probably can be responsible for an aggregate shape */
                glm::dvec3 v;
                v.x = v12.x/d;
                v.y = v12.y/d;
                v.z = v12.z/d;
                t_confinements.push_back( v );
                isModified = true;
            }
        }
    }

    size_t t_conf_num = 0;

    if (isModified){
        size_t tc_size = t_confinements.size();
        if ( tc_size > 2 ){
            for (size_t i = 0; i < tc_size; i++){
                for (size_t j = 0; j < tc_size; j++){
                    if ( i != j ){
                        double res = std::abs( glm::dot( t_confinements[i],t_confinements[j] ) );
                        if ( res < 0.05 ){
                            //confinements.push_back(t_confinements[i]);
                            t_conf_num++;
                            break;
                        }
                    }
                }
            }
        }
        else {
            t_conf_num = t_confinements.size();
        }
    }

    //conf_num = t_conf_num;
    double t_conf_degree = 1.0 - 1.0 / ( (double)(t_conf_num + 1) );

    t_confinements.clear();
    t_confinements.shrink_to_fit();

    return t_conf_degree;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Cell::getGrowthDirections( std::vector<Cell> &cells )
{
    srand (time(NULL));
    getConfinements( cells );
    if ( conf_num == 0 ){
        // rand() % 1 - 1
        glm::dvec3 nvect;
        nvect.x = (double)(rand() % 1 - 1);
        nvect.y = (double)(rand() % 1 - 1);
        nvect.z = (double)(rand() % 1 - 1);

        double rhs = glm::dot( origin,nvect );
        glm::dvec3 r1;
        glm::dvec3 r2;

        r1.x = (double)(10*(rand() % 1 - 1));
        r1.y = (double)(rand() % 1 - 1);
        r1.z = ( rhs - nvect.x * r1.x - nvect.y * r1.y )/nvect.z;
        r2.x = -r1.x;
        r2.y = -r1.y;
        r2.z = ( rhs + nvect.x * r1.x + nvect.y * r1.y )/nvect.z;
        growthdir = getNormale( origin,r1,r2 );
    }
    else if ( conf_num == 1 ){
        glm::dvec3 orthdir = getOrthVect( confinements[0] );
        growthdir = getCross( confinements[0], orthdir );
    }
    else /* if ( (conf_num > 1) && (conf_degree < growthThreshold) ) */ {
        std::vector<glm::dvec3> tmp_directions;
        std::vector<glm::dvec3> directions;
        for (size_t i = 0; i < conf_num-1; i++){
            for (size_t j = i+1; j < conf_num; j++){
                glm::dvec3 confinement1 = confinements[i];
                glm::dvec3 confinement2 = confinements[j];
                glm::dvec3 dir = getCross( confinement1,confinement2 );
                if ( glm::length(dir) != 0.0 ){
                    tmp_directions.push_back(dir);
                    //tmp_directions.push_back(-dir);
                }
            }
        }
        if ( tmp_directions.size() == 0 ){
            glm::dvec3 orthdir = getOrthVect( confinements[0] );
            tmp_directions.push_back( getCross( confinements[0], orthdir ) );
        }
        
        size_t dirNum = tmp_directions.size();
        size_t acceptedNum = 0;
        for (size_t i = 0; i < dirNum; i++){
            glm::dvec3 dir = tmp_directions[i];
            bool accepted = true;
            for (size_t j = 0; j < conf_num; j++){
                double check = glm::dot( dir,confinements[j] );
                if (check < -0.9){
                    //std::cout<<"check = "<<check<<std::endl;
                    accepted = false;
                    break;
                }
            }
            if (accepted){
                directions.push_back(dir);
                acceptedNum++;
            }
        }
        //std::cout<<"getGrowthDirections 4-4, acceptedNum = "<<acceptedNum<<std::endl; std::cout<<"getGrowthDirections -> dirNum = "<<dirNum<<std::endl; std::cout<<"conf_degree = "<<conf_degree<<std::endl;
        dirNum = directions.size();
        if (dirNum == 0){
            growthdir = glm::dvec3(0.0);
        }
        else {
            growthdir = directions[ rand() % dirNum ];
        }      
        
    }
    // else {
    //     growthdir = glm::dvec3(0.0);
    // }
}

//---------------------------------------------------------------------------------------------------------------------------------

glm::dvec3 Cell::getNewOrigin( glm::dvec3 currentorigin, glm::dvec3 growthdirection, double cellsshift )
{
    growthdirection.x = growthdirection.x * cellsshift;
    growthdirection.y = growthdirection.y * cellsshift;
    growthdirection.z = growthdirection.z * cellsshift;

    return currentorigin - growthdirection;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Cell::correct_origin()
{
        origin.x = x;
        origin.y = y;
        origin.z = z;
}

