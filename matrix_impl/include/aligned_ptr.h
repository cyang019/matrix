#ifndef _ALIGN_ALIGNED_PTR_H
#define _ALIGN_ALIGNED_PTR_H

#include <memory>
#include <type_traits>
#include <cstdlib>    // aligned_alloc
#include <complex>


namespace align {

#if defined(MSVC)
  template<typename T> using unique_ptr_aligned = std::unique_ptr<T, decltype(&_aligned_free)>;
#else
  inline void* aligned_malloc(size_t bytes, size_t alignment)
  {
      void* p1; // original block
      void* p2; // aligned block
      size_t offset = alignment - 1 + sizeof(void*);

      p1 = std::malloc((bytes + alignment - 1) + sizeof(void *));
      if(!p1)
      {
         return NULL;
      }
      p2 = (void*)(((size_t)(p1) + offset) & ~(alignment - 1));

      *((void**)p2 - 1) = p1;
      return p2;
  }
  
  inline void aligned_free(void *p)
  {
    std::free(((void**)p)[-1]);
  }

  template<typename T> using unique_ptr_aligned = std::unique_ptr<T, decltype(&aligned_free)>;
#endif

  template<typename T>
  unique_ptr_aligned<T> aligned_uptr(size_t align, size_t size)
  {
#if defined(MSVC)
      auto res = unique_ptr_aligned<T>(
          static_cast<std::remove_all_extents_t<T> *>(_aligned_malloc(align, size)), 
          &_aligned_free);
#else
      auto res = unique_ptr_aligned<T>(
          static_cast<std::remove_all_extents_t<T> *>(aligned_malloc(align, size)), 
          &aligned_free);
#endif
      return res;
  }
}


#endif
