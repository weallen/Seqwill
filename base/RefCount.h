#ifndef REFCOUNT_H_
#define REFCOUNT_H_

#include <boost/intrusive_ptr.hpp>
#include <boost/detail/atomic_count.hpp>

#include "base/Common.h"
class RefBase
{
public:
    typedef boost::intrusive_ptr<RefBase> Ptr;

    RefBase()
        : m_counter(0)
    {}

    virtual ~RefBase() 
    {}

    size_t RefCount()
    { return m_counter; }

    void AddRef() const
    {
        ++m_counter;
    }

    bool Release() const
    {
        return(--m_counter == 0);
    }

private:
    mutable boost::detail::atomic_count m_counter;

    DISALLOW_COPY_AND_ASSIGN(RefBase);
};


inline void
intrusive_ptr_add_ref(RefBase* r)
{
    r->AddRef();
}

inline void
intrusive_ptr_release(RefBase* r)
{
    r->Release();
    if (r->RefCount() == 0)
        delete r;
}
#endif
