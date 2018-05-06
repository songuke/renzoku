#include "json_parser.h"
#include "log.h"

#include <json.h>

#include <sstream>
using namespace std;

namespace Renzoku {

static int my_callback(void *userdata, int type, const char *data, uint32_t length)
{
    JsonParser *parser = (JsonParser *)userdata;
    switch (type) {
    case JSON_OBJECT_BEGIN:
        parser->begin_object();
        break;
    case JSON_ARRAY_BEGIN:
        parser->begin_array();
        break;
    case JSON_OBJECT_END:
        parser->end_object();
        break;
    case JSON_ARRAY_END:
        parser->end_array();
        break;
    case JSON_KEY:
        parser->add_key(data, length);
        break;
    case JSON_STRING:
        parser->add_value(new JsonValue(new JsonPrimitive(JsonPrimitive::JSON_PRIMITIVE_STRING, data, length)));
        break;
    case JSON_INT:
        parser->add_value(new JsonValue(new JsonPrimitive(JsonPrimitive::JSON_PRIMITIVE_INT, data, length)));
        break;
    case JSON_FLOAT:
        parser->add_value(new JsonValue(new JsonPrimitive(JsonPrimitive::JSON_PRIMITIVE_FLOAT, data, length)));
        break;
    case JSON_NULL:
        //parser->set_value_null();
        break;
    case JSON_TRUE:
        parser->add_value(new JsonValue(new JsonPrimitive(JsonPrimitive::JSON_PRIMITIVE_BOOL, "1", 1)));
        break;
    case JSON_FALSE:
        parser->add_value(new JsonValue(new JsonPrimitive(JsonPrimitive::JSON_PRIMITIVE_BOOL, "0", 1)));
        break;
    }
    return 0;
}

static int process_file(json_parser *parser, FILE *input, int *retlines, int *retcols)
{
    char buffer[4096];
    int ret = 0;
    int32_t read;
    int lines, col, i;

    lines = 1;
    col = 0;
    while (1) {
        uint32_t processed;
        read = fread(buffer, 1, 4096, input);
        if (read <= 0)
            break;
        ret = json_parser_string(parser, buffer, read, &processed);
        for (i = 0; i < processed; i++) {
            if (buffer[i] == '\n') { col = 0; lines++; }
            else col++;
        }
        if (ret)
            break;
    }
    if (retlines) *retlines = lines;
    if (retcols) *retcols = col;
    return ret;
}

JsonParser::JsonParser(JsonElementVisitor *visitor) : visitor(visitor) {}

void JsonParser::load(const string file) {
    json_parser parser;
    if (json_parser_init(&parser, NULL, my_callback, this)) {
        Log::critical("Libjson initialization went wrong.", ExitCode::EXIT_CODE_LIBRARY_FAILED_TO_LOAD);
    }

    FILE *f = fopen(file.c_str(), "rb");
    if (!f) {
        ostringstream oss;
        oss << "Cannot open " << file << ": " << strerror(errno);
        Log::critical(oss.str(), ExitCode::EXIT_CODE_FILE_FAILED_TO_LOAD);
        return;
    }

    int ret, lines, cols;
    ret = process_file(&parser, f, &lines, &cols);
    if (ret) {
        ostringstream oss;
        oss << "JSON file syntax error near line " << lines << ", column " << cols << ".";
        Log::critical(oss.str(), ExitCode::EXIT_CODE_SCENE_FILE_ERROR);
        return;
    }

    ret = json_parser_is_done(&parser);
    if (!ret) {
        Log::critical("JSON file syntax error.", ExitCode::EXIT_CODE_SCENE_FILE_ERROR);
        return;
    }

    fclose(f);
    json_parser_free(&parser);
}

void JsonParser::add_key(const char *data, int length) {
    JsonObject *obj = get_current_object();

    obj->keys.push_back(data);
}

JsonContainer *JsonParser::get_current_container() const {
    if (stack.size() > 0) {
        JsonContainer *container = stack[stack.size() - 1];
        return container;
    }
    return NULL;
}

JsonObject *JsonParser::get_current_object() const {
    if (stack.size() > 0) {
        JsonContainer *container = stack[stack.size() - 1];
        if (container->type == JsonContainer::JSON_CONTAINER_OBJECT)
            return container->data.object;
        else return NULL;
    }
    return NULL;
}

JsonArray *JsonParser::get_current_array() const {
    if (stack.size() > 0) {
        JsonContainer *container = stack[stack.size() - 1];
        if (container->type == JsonContainer::JSON_CONTAINER_ARRAY)
            return container->data.array;
        else return NULL;
    }
    return NULL;
}

string JsonParser::get_last_key() const {
    JsonObject *object = get_current_object();
    if (object && object->keys.size() > 0)
        return object->keys[object->keys.size() - 1];

    return "";
}

void JsonParser::delete_last_key() {
    JsonObject *object = get_current_object();
    if (object && object->keys.size() > 0)
        object->keys.pop_back();
}

void JsonParser::begin_object() {
    JsonObject *obj = new JsonObject();

    stack.push_back(new JsonContainer(obj));
}

void JsonParser::end_object() {
    JsonObject *obj = get_current_object();
    if (obj) {
        stack.pop_back();
        // try process this object
        string last_key = get_last_key();
        bool done = false;
        if (last_key != "")
            done = process_object(last_key, obj);

        if (!done) {
            if (stack.empty())
                process_root(obj);
            else
                add_object(obj);
        }
        else {
            delete obj;
            delete_last_key();
        }
    }
    else {
        Log::error() << "Object syntax error." << endn;
    }
}

void JsonParser::begin_array() {
    JsonArray *array = new JsonArray();
    stack.push_back(new JsonContainer(array));
}

void JsonParser::end_array() {
    JsonArray *array = get_current_array();
    if (array) {
        stack.pop_back();
        // try process this object
        string last_key = get_last_key();
        bool done = false;
        if (last_key != "")
            done = process_array(last_key, array);

        if (!done) {
            if (stack.empty())
                Log::error() << "JSON array as root object." << endn;
            else
                add_array(array);
        }
        else {
            delete array;
            delete_last_key();
        }
    }
    else {
        Log::error() << "Array syntax error." << endn;
    }
}

void JsonParser::add_value(JsonValue *value) {
    JsonContainer *container = get_current_container();
    if (container->type == JsonContainer::JSON_CONTAINER_OBJECT)
        container->data.object->values.push_back(value);
    else
        container->data.array->values.push_back(value);
}

void JsonParser::add_array(JsonArray *array) {
    add_value(new JsonValue(array));
}

void JsonParser::add_object(JsonObject *object) {
    add_value(new JsonValue(object));
}

bool JsonParser::process_object(string last_key, JsonObject *obj) {
    int stack_level = stack.size();

    return visitor->process_object(stack_level, last_key, obj);
}
bool JsonParser::process_array(string last_key, JsonArray *array) {
    return visitor->process_array(last_key, array);
}

void JsonParser::process_root(JsonObject *root) {
    visitor->process_root(root);
}

}  // end namespace