#include "aggregate_bvh.h"
#include "triangle.h"
#include "log.h"
#include "stats.h"

extern "C" {
    #include "morton.h"
}
using namespace Nvidia;

#include <queue>
using namespace std;

namespace Renzoku {

AggregateBvh::AggregateBvh(Surfaces &_surfaces) : 
    Aggregate(_surfaces), 
    sorted_triangles(NULL), sorted_morton(NULL), sorted_surface_id(NULL), sorted_boxes(NULL),
    root(NULL), nodes(NULL), num_nodes(0)
{    
    build();
    Aggregate::compute_density();
}

AggregateBvh::~AggregateBvh() {
    cleanup();    
}

void AggregateBvh::cleanup() {
    if (sorted_triangles)   delete [] sorted_triangles;
    if (sorted_morton)      delete [] sorted_morton;
    if (sorted_surface_id)  delete [] sorted_surface_id;
    if (sorted_boxes)       delete [] sorted_boxes;
    if (nodes)              delete [] nodes;
}

static void get_depth_first_order(BvhNode *nodes, int cur, int *order, int &num_nodes_so_far) {
    order[cur] = num_nodes_so_far;
    num_nodes_so_far++;

    BvhNode &node = nodes[cur];
    if (node.index_left >= 0)
        get_depth_first_order(nodes, node.index_left, order, num_nodes_so_far);
    if (node.index_right >= 0)
        get_depth_first_order(nodes, node.index_right, order, num_nodes_so_far);
}

void AggregateBvh::build() {
    // allocate triangle array
    num_triangles = 0;
    for (int i = 0; i < surfaces.size(); ++i) {
        num_triangles += surfaces[i].get_shape()->get_triangle_count();
    }
    if (num_triangles <= 0) return;

    Triangle *triangles = new Triangle[num_triangles];
    BoundingBox *boxes = new BoundingBox[num_triangles];
    int *surface_id = new int[num_triangles];
    box.reset();    // global bounding box
    int k = 0;
    for (int i = 0; i < surfaces.size(); ++i) {
        Shape *shape = surfaces[i].get_shape();
        vector<Triangle> tris;
        shape->get_triangles(tris);
        // shape->get_primitive_shapes(prims);      // no grouping  
        for (int j = 0; j < shape->get_triangle_count(); ++j) {
            triangles[k] = tris[j];            
            surface_id[k] = i;
            BoundingBox b = triangles[k].get_bounding_box();
            boxes[k] = b;
            box.merge(b);

            k++;
        }
    }

    Log::info() << "Total triangles     : " << num_triangles << endn;
    
    Log::info() << "BVH statistics:" << endn;
    Stats stats;

    Log::info() << "Building Morton code..." << endn;
    stats.tic();
    unsigned int *morton = new unsigned int[num_triangles];
    for (int i = 0; i < num_triangles; ++i) {
        morton[i] = morton_box(&boxes[i], &box);
    }
    stats.toc();
    stats.print_elapsed_milliseconds();
        
    Log::info() << "Sorting Morton array..." << endn;
    stats.tic();
    MortonEntry *entries;
    sort_morton_array(morton, num_triangles, entries);    

    // permute triangles
    sorted_triangles    = new Triangle[num_triangles];
    sorted_morton       = new unsigned int[num_triangles];
    sorted_surface_id   = new int[num_triangles];
    sorted_boxes        = new BoundingBox[num_triangles];
    for (int i = 0; i < num_triangles; ++i) {                
        sorted_morton[i]        = entries[i].code;
        int original_idx        = entries[i].index;
        sorted_triangles[i]     = triangles[original_idx];        
        sorted_surface_id[i]    = surface_id[original_idx];
        sorted_boxes[i]         = boxes[original_idx];
    }
    delete [] triangles;
    delete [] morton;
    delete [] surface_id;
    delete [] boxes;
    delete [] entries;
    stats.toc();
    stats.print_elapsed_milliseconds();

    Log::info() << "Building levels..." << endn;
    stats.tic();
    // build all levels    
    int num_estimated_levels = (int)ceil(log(num_triangles) / log(2.0f));    // is it a good bound?
    int capacity = (int)pow(2, num_estimated_levels + 1);
    nodes = new BvhNode[capacity];
    
    // the queue
    int first = 0,                      // current to-pop node
        last = 1,                       // next empty node
        next_last = 1;                  // next level

    num_nodes = 1;
    root = nodes;
    root->start = 0;
    root->end = num_triangles - 1;      // inclusive
    
    typedef pair<int, int> Level;       // record level info for bottom-up box building
    vector<Level> levels;
    levels.push_back(Level(0, 1));
    int level = 0;
    while (true) {        
        BvhNode* node = &nodes[first]; 
        node->level = level;

        int start = node->start;
        int end = node->end;
        if (start == end) {
            node->index_left = -1;
            node->index_right = -1;

        } else {

            int split = Nvidia::findSplit(sorted_morton, start, end);
            
            node->index_left    = next_last;
            node->index_right   = next_last + 1;
            next_last += 2;
            num_nodes += 2;            

            BvhNode *left  = &nodes[node->index_left];
            BvhNode *right = &nodes[node->index_right];
            left->start     = start;
            left->end       = split;
            right->start    = split + 1;
            right->end      = end; 
        }

        first++;
        if (first == last) {                // level end
            if (next_last == last) {
                break;                      // no nodes queued, done
            } else {                
                last = next_last;           // go to next level
                levels.push_back(Level(first, last));
                //Log::info() << "New level : " << first << " " << last << endn;
                //Log::info() << "Num nodes: " << num_nodes << endn;
                level++;
            }
        }
    }
    Log::info() << "BVH nodes/capacity : " << num_nodes << "/" << capacity << endn;
    Log::info() << "BVH levels         : " << levels.size() << endn;
    num_levels = levels.size();
    stats.toc();
    stats.print_elapsed_milliseconds();

    Log::info() << "Linking levels..." << endn;
    stats.tic();
    for (int i = levels.size() - 1; i >= 0; --i) {
        int first = levels[i].first;
        int last = levels[i].second;

        for (int j = first; j < last; ++j) {
            BvhNode *node = &nodes[j];

            if (node->index_left < 0) {
                if (node->index_right < 0) {
                    node->box = sorted_boxes[node->start];                    
                } else {
                    node->box = nodes[node->index_right].box;                    
                }
            } else {
                if (node->index_right < 0) {
                    node->box = nodes[node->index_left].box;                    
                } else {
                    node->box = nodes[node->index_left].box;
                    node->box.merge(nodes[node->index_right].box);                    
                }
            }
        }
    }
    stats.toc();
    stats.print_elapsed_milliseconds();

    //Log::debug() << "BVH bounding box : min = " << root->box.v_min << ", max = " << root->box.v_max << endn;

    /**
     * Re-organize nodes into DFS order.
     */
    // NOTE: this is slightly slower than the BFS order.
    /*
    int *dfs_order = new int[num_nodes];    // map an old node to a new node index
    int num_nodes_so_far = 0;
    get_depth_first_order(nodes, 0, dfs_order, num_nodes_so_far);
    BvhNode *dfs_nodes = new BvhNode[num_nodes];
    for (int i = 0; i < num_nodes; ++i) {
        int j = dfs_order[i];

        dfs_nodes[j] = nodes[i];
        if (nodes[i].index_left >= 0)
            dfs_nodes[j].index_left = dfs_order[nodes[i].index_left];
        else
            dfs_nodes[j].index_left = -1;

        if (nodes[i].index_right >= 0)
            dfs_nodes[j].index_right = dfs_order[nodes[i].index_right];
        else
            dfs_nodes[j].index_right = -1;
    }
    delete [] nodes;
    nodes = dfs_nodes;
    root = nodes;
    */
}


struct Stack {
    int buffer[64];
    int capacity;
    int top;

    Stack() {        
        this->capacity = 64;
        top = -1;
    }
    ~Stack() {
        
    }

    void push(int val) {
        if (top < capacity - 1) {
            ++top;
            buffer[top] = val;
        } else {
            Log::info() << "Stack out of capacity." << endn;
        }
    }

    int pop() {
        if (top >= 0) {
            int val = buffer[top];
            top--;
            return val;
        }
        return -1;
    }

    bool is_empty() {
        return top < 0;
    }

    int size() {
        return top + 1;
    }
};

struct CircularQueue {
    int *buffer;
    int capacity;
    int first, last;
    int num_items;

    CircularQueue(int capacity) {
        this->buffer = new int[capacity];
        this->capacity = capacity;
        first = 0;
        last = 0;
        num_items = 0;
    }
    ~CircularQueue() {
        delete [] buffer;
    }

    void push(int val) {
        if (num_items < capacity) {
            buffer[last] = val;
            last = (last + 1) % capacity;
            num_items++;
        } else {
            Log::info() << "Circular queue out of capacity." << endn;
        }
    }

    int pop() {
        if (num_items > 0) {
            int val = buffer[first];
            first = (first + 1) % capacity;
            num_items--;
            return val;
        }
        return -1;
    }

    bool is_empty() {
        return num_items == 0;
    }
};

/*
static bool test_hit(const Ray &r, BvhNode *nodes, Triangle *sorted_triangles, int cur, Float tmin, Float tmax, Float time, GeometryHit &record, int &tri_index) {
    record.hit = false;
    BvhNode &node = nodes[cur];    
    if (! node.box.hit(r, tmin, tmax)) return false;
    
    if (node.index_left < 0 && node.index_right < 0) {                        
        
        Triangle &tri = sorted_triangles[node.start];        
        tri.hit(r, tmin, tmax, time, record);
        if (record.hit) tri_index = node.start;
        return record.hit;

    }  else {
                
        if (node.index_left >= 0) 
            if (test_hit(r, nodes, sorted_triangles, node.index_left, tmin, tmax, time, record, tri_index)) {
                tmax = record.t;    // cut tmax down early, improve speed by 0.1 samples/sec.
            }

        if (node.index_right >= 0) {
            GeometryHit rec_right;
            int tri_right;
            if (test_hit(r, nodes, sorted_triangles, node.index_right, tmin, tmax, time, rec_right, tri_right)) {            
                record = rec_right;
                tri_index = tri_right;
            }
        }        
    }    
    return record.hit;
}

bool AggregateBvh::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    record.hit = false;
    if (num_nodes == 0) return false;

    int triangle_index;
    GeometryHit gh;
    test_hit(r, nodes, sorted_triangles, 0, tmin, tmax, time, gh, triangle_index);
    if (gh.hit) {
        Triangle *t = &sorted_triangles[triangle_index];
        t->fill_hit_record(r, gh, record);
        Surface *s = surfaces[sorted_surface_id[triangle_index]];
        s->fill_hit_record(r, record);
    }
    return record.hit;
}*/

/*
static bool test_hit2(const Ray &r, BvhNode *nodes, Triangle *sorted_triangles, int cur, Float tmin, Float tmax, Float time) {
    BvhNode &node = nodes[cur];    
    if (! node.box.hit(r, tmin, tmax)) return false;
    
    if (node.index_left < 0 && node.index_right < 0) {                        

        Triangle &tri = sorted_triangles[node.start];            
        return tri.hit(r, tmin, tmax, time);

    }  else {
                
        if (node.index_left >= 0) 
            if (test_hit2(r, nodes, sorted_triangles, node.index_left, tmin, tmax, time)) return true;

        if (node.index_right >= 0) {            
            return test_hit2(r, nodes, sorted_triangles, node.index_right, tmin, tmax, time);
        }
        
    }
    return false;
}

bool AggregateBvh::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    if (num_nodes == 0) return false;
    return test_hit2(r, nodes, sorted_triangles, 0, tmin, tmax, time);    
}
*/

bool AggregateBvh::hit(const Ray &r, Float tmin, Float tmax, Float time, HitRecord &record) const {
    if (! root) return false;
    if (! root->box.hit(r, tmin, tmax)) return false;

    //CircularQueue q(1024);
    Stack q;
    q.push(0);

    //Log::info() << "---------------" << endn;
    bool hit = false;
    int triangle_index;
    GeometryHit best_gh;
    while (! q.is_empty()) {
        int node_index = q.pop();
        BvhNode &node = nodes[node_index];
        //Log::info() << "Check node : " << node_index << " " << node.box.v_min << " / " << node.box.v_max << endn;
        if (node.index_left < 0 && node.index_right < 0) {                        

            Triangle &tri = sorted_triangles[node.start];            
            GeometryHit gh;
            if (tri.hit(r, tmin, tmax, time, gh)) {
                tmax = gh.t;                // update the best tmax so far to prune quickly
                best_gh = gh;
                hit = true;
                triangle_index = node.start;

                //Log::info() << "Triange : " << tmax << endn;
            }

        }  else {

            if (node.index_left >= 0 && nodes[node.index_left].box.hit(r, tmin, tmax)) 
                q.push(node.index_left);
            if (node.index_right >= 0 && nodes[node.index_right].box.hit(r, tmin, tmax))
                q.push(node.index_right);

            //Log::info() << "Stack : " << q.size() << endn;
        }
    }
    if (hit) {
        Triangle &t = sorted_triangles[triangle_index];
        t.fill_hit_record(r, best_gh, record);
        const Surface &s = surfaces[sorted_surface_id[triangle_index]];
        s.fill_hit_record(r, record);
    }
    return hit;
}



bool AggregateBvh::hit(const Ray &r, Float tmin, Float tmax, Float time) const {
    if (! root) return false;
    if (! root->box.hit(r, tmin, tmax)) return false;

    //CircularQueue q(1024);
    Stack q;
    q.push(0);
        
    while (! q.is_empty()) {
        BvhNode &node = nodes[q.pop()];        

        if (node.index_left < 0 && node.index_right < 0) {                        

            Triangle &tri = sorted_triangles[node.start];            
            if (tri.hit(r, tmin, tmax, time)) return true;

        }  else {
            
            if (node.index_left >= 0 && nodes[node.index_left].box.hit(r, tmin, tmax))
                q.push(node.index_left);
            if (node.index_right >= 0 && nodes[node.index_right].box.hit(r, tmin, tmax))
                q.push(node.index_right);

        }
    }
    return false;
}

} // end namespace Renzoku
