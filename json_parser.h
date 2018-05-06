#ifndef _JSON_PARSER_H_
#define _JSON_PARSER_H_

#include "math3.h"

#include <vector>
#include <string>
#include <cstring>
#include <cerrno>

using namespace std;

namespace Renzoku {

/**
* There are two styles to write a JSON parser:
* - DOM style: read the entire JSON file into a tree. Traverse the tree from the root in any orders to process the data.
*   DOM is easy to implement, provide nice and clear code, but can cost memory as the entire tree needs to be stored.
*
* - SAX style: event-based; every time an opening/close tag is detected, a callback or delegate function is called.
*   In the parser we maintain a stack of JsonObjects. The top of the stack indicates the current object which is receiving
*   data from the parser.
*   When a close tag is encountered, the top stack object is completed, and can be processed further, and/or attached to
*   the next object in the stack.
*   If we don't count the processed data, the stack depth is equal to the JSON object depth.
*
* This parser experiments with the SAX style.
*/

struct JsonException : std::exception {
    std::string s;
    JsonException(std::string ss) : s(ss) {}
    ~JsonException() throw () {}
    const char* what() const throw() { return s.c_str(); }
};

struct JsonPrimitive {
    enum JsonPrimitiveType {
        JSON_PRIMITIVE_STRING,
        JSON_PRIMITIVE_INT,
        JSON_PRIMITIVE_FLOAT,
        JSON_PRIMITIVE_BOOL
    };

    JsonPrimitive(const char *data, int length) : length(length), type(JSON_PRIMITIVE_STRING) {
        this->data = new char[length + 1];
        strcpy(this->data, data);
    }
    JsonPrimitive(JsonPrimitiveType type, const char *data, int length) : length(length), type(type) {
        this->data = new char[length + 1];
        strcpy(this->data, data);
    }
    ~JsonPrimitive() {
        if (data)
            delete[] data;
    }

    char *data;
    int length;
    int type;

    bool is_int() const {
        return type == JSON_PRIMITIVE_INT;
    }

    bool is_float() const {
        return type == JSON_PRIMITIVE_FLOAT;
    }

    bool is_string() const {
        return type == JSON_PRIMITIVE_STRING;
    }

    int to_int() {
        if (type == JSON_PRIMITIVE_INT) {
            return atoi(data);
        }
        throw JsonException("Invalid type cast.");
    }

    float to_float() {
        if (type == JSON_PRIMITIVE_FLOAT) {
            return (float)atof(data);
        }
        throw JsonException("Invalid type cast.");
    }

    double to_double() {
        if (type == JSON_PRIMITIVE_FLOAT) {
            return atof(data);
        }
        throw JsonException("Invalid type cast.");
    }

    bool to_bool() {
        if (type == JSON_PRIMITIVE_BOOL) {
            if (data[0] == '0')
                return false;
            else
                return true;
        }
        throw JsonException("Invalid type cast.");
    }

    string to_string() {
        if (type == JSON_PRIMITIVE_STRING) {
            return data;
        }
        throw JsonException("Invalid type cast.");
    }
};


struct JsonObject;
struct JsonArray;

struct JsonValue {
    enum JsonValueType {
        JSON_VALUE_PRIMITIVE,
        JSON_VALUE_OBJECT,
        JSON_VALUE_ARRAY
    };
    int type;

    union {
        JsonPrimitive *primitive;
        JsonObject *object;
        JsonArray *array;
    } data;

    JsonValue(JsonObject *object) {
        type = JSON_VALUE_OBJECT;
        data.object = object;
    }

    JsonValue(JsonPrimitive *primitive) {
        type = JSON_VALUE_PRIMITIVE;
        data.primitive = primitive;
    }

    JsonValue(JsonArray *array) {
        type = JSON_VALUE_ARRAY;
        data.array = array;
    }

    JsonArray *to_array() {
        if (type == JSON_VALUE_ARRAY)
            return data.array;

        throw JsonException("Invalid type cast.");
    }

    JsonObject *to_object() {
        if (type == JSON_VALUE_OBJECT)
            return data.object;

        throw JsonException("Invalid type cast.");
    }

    bool is_string() {
        if (type == JSON_VALUE_PRIMITIVE && data.primitive->type == JsonPrimitive::JSON_PRIMITIVE_STRING) {
            return true;
        }
        return false;
    }

    bool is_primitive() {
        if (type == JSON_VALUE_PRIMITIVE)
            return true;
        return false;
    }

    bool is_object() {
        if (type == JSON_VALUE_OBJECT)
            return true;
        return false;
    }

    bool is_array() {
        if (type == JSON_VALUE_ARRAY)
            return true;
        return false;
    }

    string to_string() {
        if (type == JSON_VALUE_PRIMITIVE)
            return data.primitive->to_string();

        throw JsonException("Invalid type cast.");
    }

    int to_int() {
        if (type == JSON_VALUE_PRIMITIVE)
            return data.primitive->to_int();

        throw JsonException("Invalid type cast.");
    }

    float to_numeric() {
        if (type == JSON_VALUE_PRIMITIVE) {
            if (data.primitive->is_int())
                return (float)data.primitive->to_int();
            if (data.primitive->is_float())
                return (float)data.primitive->to_float();

            throw JsonException("Invalid type cast.");
        }
        throw JsonException("Invalid type cast.");
    }

    bool to_bool() {
        if (type == JSON_VALUE_PRIMITIVE)
            return data.primitive->to_bool();

        throw JsonException("Invalid type cast.");
    }
};

struct JsonArray {
    vector<JsonValue *> values;

    Vec2 to_vec2() const {
        return Vec2(values[0]->to_numeric(),
            values[1]->to_numeric());
    }

    Vec3 to_vec3() const {
        return Vec3(values[0]->to_numeric(),
            values[1]->to_numeric(),
            values[2]->to_numeric());
    }

    Rgb to_rgb() const {
        return Rgb(values[0]->to_numeric(),
            values[1]->to_numeric(),
            values[2]->to_numeric());
    }

    Size2 to_size2() const {
        return Size2(values[0]->to_numeric(),
            values[1]->to_numeric());
    }
};

struct JsonObject {
    vector<string> keys;
    vector<JsonValue *> values;

    bool is_empty() const {
        return (keys.size() == 0 && values.size() == 0);
    }

    bool is_empty_keys() const {
        return keys.size() == 0;
    }

    bool is_empty_values() const {
        return values.size() == 0;
    }

    JsonValue *get_value(string key) {
        for (int i = 0; i < keys.size(); ++i)               // only a few pairs, use linear search
            if (key.compare(keys[i]) == 0)
                return values[i];
        return NULL;
    }

    //
    // Some utility functions to retrieve values.
    //

    /**
    * Use this function if the key is compulsory.
    */
    JsonValue *get_value_strict(string key) {
        JsonValue *val = get_value(key);
        if (!val) {
            throw JsonException("Key " + key + " not found.");
        }
        return val;
    }

    bool has_key(string key) {
        for (int i = 0; i < keys.size(); ++i)               // only a few pairs, use linear search
            if (key.compare(keys[i]) == 0)
                return true;
        return false;
    }

    int get_int(string key) {
        return get_value_strict(key)->to_int();
    }

    bool get_bool(string key) {
        return get_value_strict(key)->to_bool();
    }

    float get_numeric(string key) {
        return get_value_strict(key)->to_numeric();
    }

    string get_string(string key) {
        return get_value_strict(key)->to_string();
    }

    Vec2 get_vec2(string key) {
        return get_value_strict(key)->to_array()->to_vec2();
    }

    Vec3 get_vec3(string key) {
        return get_value_strict(key)->to_array()->to_vec3();
    }

    Rgb get_rgb(string key) {
        return get_value_strict(key)->to_array()->to_rgb();
    }

    Size2 get_size2(string key) {
        return get_value_strict(key)->to_array()->to_size2();
    }
};

struct JsonContainer {
    enum JsonContainerType {
        JSON_CONTAINER_OBJECT,
        JSON_CONTAINER_ARRAY
    };

    int type;
    union {
        JsonObject *object;
        JsonArray *array;
    } data;

    JsonContainer(JsonObject *object) {
        type = JSON_CONTAINER_OBJECT;
        data.object = object;
    }

    JsonContainer(JsonArray *array) {
        type = JSON_CONTAINER_ARRAY;
        data.array = array;
    }
};

class JsonElementVisitor {
public:
    virtual bool process_object(int stack_level, string last_key, JsonObject *obj) = 0;
    virtual bool process_array(string last_key, JsonArray *array) = 0;
    virtual void process_root(JsonObject *root) = 0;
};

/**
* This class is a bridge between the low-level events in libjson to retrieve the JSON data
* and a "visitor" or "importer" class which is responsible to process the JSON array and objects found.
* e.g., load the scene from such data.
*/
class JsonParser {
public:
    JsonParser(JsonElementVisitor *visitor);
    void load(const string file);

    void begin_object();
    void end_object();

    void begin_array();
    void end_array();

    void add_key(const char *data, int length);
    void add_value(JsonValue *value);

private:
    JsonContainer *get_current_container() const;
    JsonObject *get_current_object() const;
    JsonArray *get_current_array() const;

    string get_last_key() const;
    void delete_last_key();

    /**
    * Try to process an object or array as soon as it is ready.
    */
    bool process_object(string last_key, JsonObject *obj);
    bool process_array(string last_key, JsonArray *array);

    /**
    * The tree (similar to DOM) is eventually created when the end object is the root object.
    *
    * This function can work on the unprocessed objects left in the tree.
    */
    void process_root(JsonObject *root);

    /**
    * Unprocessed objects and arrays are handed back to the parent object for processing in later attempts.
    */
    void add_object(JsonObject *obj);
    void add_array(JsonArray *array);

private:
    vector<JsonContainer *> stack;

    JsonElementVisitor *visitor;
};

}  // end namespace

#endif