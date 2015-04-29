/**

* Author : Nguyen Minh Thang
* 
 *
 */

#include <tulip/Algorithm.h>
#include <tulip/TulipPluginHeaders.h>
#include <tulip/Delaunay.h>
#include <tulip/LayoutProperty.h>
#include <queue>
#include <math.h>
#include <cmath>
#include <algorithm>
 using namespace std;

 class PriorityEdge
 {
 public:
    PriorityEdge(double _priority,pair<unsigned int,unsigned int> _edge){
        priority = _priority;
        edge = _edge;
    }
    double getPriority(){
        return priority;
    }
    pair<tlp::node, tlp::node> getEdge(){
        return edge;
    }
private:
    double priority;
    pair<tlp::node, tlp::node> edge;
};
struct ComparePriorityEdge{
public:
    bool operator()(PriorityEdge& prt1, PriorityEdge& prt2){
        if(prt1.getPriority() >= prt2.getPriority()) return true;
        return false;
    }    
} comparePriorityEdge;
std::vector<tlp::node> iteratorToVector(tlp::Iterator<tlp::node> *iterator){
  std::vector<tlp::node> v;
  while(iterator->hasNext()){
    v.push_back(iterator->next());
}
return v;
}
void buildGraph(tlp::Graph *graph) {
  tlp::node n1=graph->addNode();
  tlp::node n2=graph->addNode();
  tlp::node n3=graph->addNode();
  tlp::node n4=graph->addNode();

  graph->addEdge(n1,n2);
  graph->addEdge(n1,n3);
  graph->addEdge(n1,n4);
  graph->addEdge(n2,n3);
  graph->addEdge(n2,n4);
}
std::vector<pair<tlp::node, tlp::node> > FindSetOfEdgesOfQuadrilateral(tlp::node p1, tlp::node p2, tlp::Graph *graph){
    std::vector<pair<tlp::node, tlp::node> > setOfEdgesOfQuadrilateral;
    tlp::Iterator<tlp::node> *neighborhoodP1 = graph->getInOutNodes(p1);
    tlp::Iterator<tlp::node> *neighborhoodP2 = graph->getInOutNodes(p2);
    std::vector<tlp::node> vP1 = iteratorToVector(neighborhoodP1);
    std::vector<tlp::node> vP2 = iteratorToVector(neighborhoodP2);
    for (size_t i = 0; i < vP1.size(); ++i)
    {
      for (size_t j = 0; j < vP2.size(); ++j)
      {
        if(vP1[i]==vP2[j]){
            setOfEdgesOfQuadrilateral.push_back(make_pair(p1,vP1[i]));
            setOfEdgesOfQuadrilateral.push_back(make_pair(p2,vP1[i]));  
        }
    }
}
return setOfEdgesOfQuadrilateral;
}
std::vector<tlp::node> FindSetPointOfQuadrilateral(tlp::node p1, tlp::node p2, tlp::Graph *graph){
    std::vector<tlp::node> setPointOfQuadrilateral;
    setPointOfQuadrilateral.push_back(p1);
    setPointOfQuadrilateral.push_back(p2);
    tlp::Iterator<tlp::node> *neighborhoodP1 = graph->getInOutNodes(p1);
    tlp::Iterator<tlp::node> *neighborhoodP2 = graph->getInOutNodes(p2);
    std::vector<tlp::node> vP1 = iteratorToVector(neighborhoodP1);
    std::vector<tlp::node> vP2 = iteratorToVector(neighborhoodP2);
    for (size_t i = 0; i < vP1.size(); ++i)
    {
      for (size_t j = 0; j < vP2.size(); ++j)
      {
        if(vP1[i]==vP2[j]){
          setPointOfQuadrilateral.push_back(vP1[i]);
      }
  }
}
return setPointOfQuadrilateral;
}
double FindRadius(double a,double b,double c){
    double radius = (a*b*c)/sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));
    return radius;
}
pair<tlp::node, tlp::node> FindFlip(tlp::node p1, tlp::node p2, tlp::Graph *graph){
    std::vector<tlp::node> flipNode; 
    tlp::Iterator<tlp::node> *neighborhoodP1 = graph->getInOutNodes(p1);
    tlp::Iterator<tlp::node> *neighborhoodP2 = graph->getInOutNodes(p2);
    std::vector<tlp::node> vP1 = iteratorToVector(neighborhoodP1);
    std::vector<tlp::node> vP2 = iteratorToVector(neighborhoodP2);
    for (size_t i = 0; i < vP1.size(); ++i)
    {
      for (size_t j = 0; j < vP2.size(); ++j)
      {
        if(vP1[i]==vP2[j]){
          flipNode.push_back(vP1[i]);
      }
  }
}
return make_pair(flipNode[0],flipNode[1]);
}
double FindEventPoint(tlp::LayoutProperty *layout,std::vector<tlp::node> setPointOfQuadrilateral,double error){
    double t = 1;
    tlp::Graph *subGraph = tlp::newGraph();
    buildGraph(subGraph); 
    std::vector<tlp::node> nodesSG;
    nodesSG.reserve(subGraph->numberOfNodes());
    tlp::LayoutProperty *layoutSubGraph = subGraph->getProperty<tlp::LayoutProperty>("viewLayout");
    tlp::node n;
    forEach(n, subGraph->getNodes()) {
        nodesSG.push_back(n);
    }
    int i = 0;
    forEach(n, subGraph->getNodes()) {
        layoutSubGraph->setNodeValue(n,tlp::Coord(round(layout->getNodeValue(setPointOfQuadrilateral[i]).getX()),layout->getNodeValue(setPointOfQuadrilateral[i]).getY(),0));
        i++;
    }
    double edgeLength = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[1]));
    double a = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[2]));
    double b = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[1],nodesSG[2]));
    double c = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[3]));
    double d = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[1],nodesSG[3]));
    double r1 = FindRadius(a,b,edgeLength);
    double r2 = FindRadius(c,d,edgeLength);
    if(edgeLength>r2 && r2>=r1){
      while(r2>r1){
          t = t+ error;
          layoutSubGraph->scale(tlp::Coord(static_cast<float>(t),static_cast<float>(1),static_cast<float>(1)));
          edgeLength = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[1]));
          a = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[2]));
          b = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[1],nodesSG[2]));
          c = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[3]));
          d = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[1],nodesSG[3]));
          r1 = FindRadius(a,b,edgeLength);
          r2 = FindRadius(c,d,edgeLength);
      }
  }else if (edgeLength>r1 && r1>r2){
      while(r2<r1){
          t = t+ error;
          layoutSubGraph->scale(tlp::Coord(static_cast<float>(t),static_cast<float>(1),static_cast<float>(1)));
          edgeLength = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[1]));
          a = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[2]));
          b = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[1],nodesSG[2]));
          c = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[0],nodesSG[3]));
          d = layoutSubGraph->edgeLength(subGraph->existEdge(nodesSG[1],nodesSG[3]));
          r1 = FindRadius(a,b,edgeLength);
          r2 = FindRadius(c,d,edgeLength);
      }
  }
  delete subGraph;
  return t;
}
bool findAspectRatio(tlp::Graph *graph,double error) {
    vector<tlp::node> nodes;
    nodes.reserve(graph->numberOfNodes());
    vector<tlp::Coord> points;
    points.reserve(graph->numberOfNodes());
    tlp::node n;
    tlp::LayoutProperty *layout = graph->getProperty<tlp::LayoutProperty>("viewLayout");
    forEach(n, graph->getNodes()) {
        nodes.push_back(n);
        points.push_back(layout->getNodeValue(n));
    }
    vector<pair<unsigned int, unsigned int> > edges;
    vector<vector<unsigned int> > simplices;
    bool ret = tlp::delaunayTriangulation(points, edges, simplices);
    if (ret) {

        graph->addCloneSubGraph("Original graph");
        
        tlp::Graph *delaunaySg = graph->addCloneSubGraph("Delaunay");
        delaunaySg->delEdges(graph->getEdges());
        
        for (size_t i = 0 ; i < edges.size() ; ++i) {
            delaunaySg->addEdge(nodes[edges[i].first], nodes[edges[i].second]);
        }

        double s = 1;
        std::vector<PriorityEdge> Q;
        for (size_t i = 0 ; i < edges.size() ; ++i) {
            std::vector<tlp::node> setPointOfQuadrilateral = FindSetPointOfQuadrilateral(nodes[edges[i].first],nodes[edges[i].second],delaunaySg);
            if(setPointOfQuadrilateral.size()==4){
                double t = FindEventPoint(layout,setPointOfQuadrilateral,error);
                if(t>s){
                    PriorityEdge priorityEdge(t,make_pair(nodes[edges[i].first],nodes[edges[i].second]));
                    Q.push_back(priorityEdge);
                }      
            }
        }
        std::sort(Q.begin(),Q.end(),comparePriorityEdge);
        while(Q.empty()==false){
          PriorityEdge ptr = Q.back();
          Q.pop_back();
          if(ptr.getPriority()>=s){
            s = ptr.getPriority();
          }
          pair<tlp::node, tlp::node> e = ptr.getEdge();
          pair<tlp::node, tlp::node> flip = FindFlip(e.first,e.second,delaunaySg);
          delaunaySg->delEdge(delaunaySg->existEdge(e.second,e.first,false));
          delaunaySg->addEdge(flip.first,flip.second);
          std::vector<pair<tlp::node, tlp::node> > quaEdges = FindSetOfEdgesOfQuadrilateral(flip.first,flip.second,delaunaySg);
          for (size_t i = 0 ; i < quaEdges.size() ; ++i) {
            tlp::edge quaEdge = delaunaySg->existEdge(quaEdges[i].first,quaEdges[i].second,false);
            int index = -1;
            for (size_t j = 0; j < Q.size(); ++j)
            {
              tlp::edge qEdge = delaunaySg->existEdge(Q[j].getEdge().first,Q[j].getEdge().second,false);
              if(quaEdge==qEdge){
                std::vector<tlp::node> setPointOfQuadrilateral = FindSetPointOfQuadrilateral(quaEdges[i].first,quaEdges[i].second,delaunaySg);
                double t = FindEventPoint(layout,setPointOfQuadrilateral,error);
                    if(t>1){
                        PriorityEdge priorityEdge(t,make_pair(quaEdges[i].first,quaEdges[i].second));
                        Q[j] = priorityEdge;
                    }else{
                      index = j;
                    }
                    break;   
              }
            }
            if(index != -1){
              Q.erase(Q.begin()+index);
            }
          }
        }
        layout->scale(tlp::Coord(static_cast<float>(s),static_cast<float>(1),static_cast<float>(1)));
        cout<<"delaunaySg"<<endl;
        cout<<delaunaySg<<flush;       
        graph->delEdges(graph->getEdges());
    }
    return ret;
}

class AspectRatio : public tlp::Algorithm {

    public :
    
    AspectRatio(tlp::PluginContext *context) : Algorithm(context) {
        addInParameter<double>("error", "This is a error parameter, it must be small than 0.05", "0.05");
    }
    
    PLUGININFORMATION("Aspect Ratio of a Scatter Plot","Quaywin","","","1.0","Triangulation")
    
    bool run() {
        tlp::Observable::holdObservers();
        double error = 0.05;
        if (dataSet)
            dataSet->get("error", error);
        bool ret = findAspectRatio(graph,error);
        
        tlp::Observable::unholdObservers();
        
        return ret;
    }
    
};

PLUGIN(AspectRatio)