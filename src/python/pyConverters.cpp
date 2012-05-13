#include "pyConverters.h"

void exportConverters()
{
    Dvector_from_PyObject();
}

Dvector_from_PyObject::Dvector_from_PyObject()
{
    converter::registry::push_back(&convertible, &construct, type_id<dvector>());
}

void* Dvector_from_PyObject::convertible(PyObject* obj_ptr)
{
    if (PySequence_Check(obj_ptr)) {
        return obj_ptr;
    } else {
        return NULL;
    }
}

void Dvector_from_PyObject::construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data)
{
    typedef converter::rvalue_from_python_storage<dvector> storage_t;
    storage_t* the_storage = reinterpret_cast<storage_t*>(data);
    void* memory_chunk = the_storage->storage.bytes;
    object obj( handle<>( borrowed(obj_ptr)));
    dvector* v = new (memory_chunk) dvector(len(obj));
    for (size_t i=0; i<v->size(); i++) {
        (*v)[i] = extract<double>(obj[i]);
    }
    data->convertible = memory_chunk;
}
