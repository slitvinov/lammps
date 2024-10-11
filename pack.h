struct pack_plan_3d {
  int nfast;
  int nmid;
  int nslow;
  int nstride_line;
  int nstride_plane;
  int nqty;
};
#if !defined(FFT_PACK_POINTER) && !defined(FFT_PACK_MEMCPY)
#define FFT_PACK_ARRAY 
#endif
#ifndef PACK_DATA
#define PACK_DATA double
#endif
#ifdef FFT_PACK_ARRAY
static void pack_3d(PACK_DATA *data, PACK_DATA *buf, struct pack_plan_3d *plan)
{
  int in, out, fast, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  in = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      out = plane + mid * nstride_line;
      for (fast = 0; fast < nfast; fast++) buf[in++] = data[out++];
    }
  }
}
static void unpack_3d(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  int in, out, fast, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + mid * nstride_line;
      for (fast = 0; fast < nfast; fast++) data[in++] = buf[out++];
    }
  }
}
static void unpack_3d_permute1_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  int in, out, fast, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + mid;
      for (fast = 0; fast < nfast; fast++, in += nstride_plane) data[in] = buf[out++];
    }
  }
}
static void unpack_3d_permute1_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  int in, out, fast, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + 2 * mid;
      for (fast = 0; fast < nfast; fast++, in += nstride_plane) {
        data[in] = buf[out++];
        data[in + 1] = buf[out++];
      }
    }
  }
}
static void unpack_3d_permute1_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  int in, out, iqty, instart, fast, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane, nqty;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;
  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      instart = plane + nqty * mid;
      for (fast = 0; fast < nfast; fast++, instart += nstride_plane) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) data[in++] = buf[out++];
      }
    }
  }
}
static void unpack_3d_permute2_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  int in, out, fast, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      in = slow + mid * nstride_plane;
      for (fast = 0; fast < nfast; fast++, in += nstride_line) data[in] = buf[out++];
    }
  }
}
static void unpack_3d_permute2_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  int in, out, fast, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      in = 2 * slow + mid * nstride_plane;
      for (fast = 0; fast < nfast; fast++, in += nstride_line) {
        data[in] = buf[out++];
        data[in + 1] = buf[out++];
      }
    }
  }
}
static void unpack_3d_permute2_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  int in, out, iqty, instart, fast, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, nqty;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;
  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      instart = nqty * slow + mid * nstride_plane;
      for (fast = 0; fast < nfast; fast++, instart += nstride_line) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) data[in++] = buf[out++];
      }
    }
  }
}
#endif
#ifdef FFT_PACK_POINTER
static void pack_3d(PACK_DATA *data, PACK_DATA *buf, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  in = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane + mid * nstride_line]);
      end = begin + nfast;
      for (out = begin; out < end; out++) *(in++) = *out;
    }
  }
}
static void unpack_3d(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane + mid * nstride_line]);
      end = begin + nfast;
      for (in = begin; in < end; in++) *in = *(out++);
    }
  }
}
static void unpack_3d_permute1_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane + mid]);
      end = begin + nfast * nstride_plane;
      for (in = begin; in < end; in += nstride_plane) *in = *(out++);
    }
  }
}
static void unpack_3d_permute1_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane + 2 * mid]);
      end = begin + nfast * nstride_plane;
      for (in = begin; in < end; in += nstride_plane) {
        *in = *(out++);
        *(in + 1) = *(out++);
      }
    }
  }
}
static void unpack_3d_permute1_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *instart, *begin, *end;
  int iqty, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane, nqty;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane + nqty * mid]);
      end = begin + nfast * nstride_plane;
      for (instart = begin; instart < end; instart += nstride_plane) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
      }
    }
  }
}
static void unpack_3d_permute2_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[slow + mid * nstride_plane]);
      end = begin + nfast * nstride_line;
      for (in = begin; in < end; in += nstride_line) *in = *(out++);
    }
  }
}
static void unpack_3d_permute2_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[2 * slow + mid * nstride_plane]);
      end = begin + nfast * nstride_line;
      for (in = begin; in < end; in += nstride_line) {
        *in = *(out++);
        *(in + 1) = *(out++);
      }
    }
  }
}
static void unpack_3d_permute2_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *instart, *begin, *end;
  int iqty, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, nqty;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[nqty * slow + mid * nstride_plane]);
      end = begin + nfast * nstride_line;
      for (instart = begin; instart < end; instart += nstride_line) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
      }
    }
  }
}
#endif
#ifdef FFT_PACK_MEMCPY
static void pack_3d(PACK_DATA *data, PACK_DATA *buf, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out;
  int mid, slow, size;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane, upto;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  size = nfast * sizeof(PACK_DATA);
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_plane;
    upto = slow * nmid * nfast;
    for (mid = 0; mid < nmid; mid++) {
      in = &(buf[upto + mid * nfast]);
      out = &(data[plane + mid * nstride_line]);
      memcpy(in, out, size);
    }
  }
}
static void unpack_3d(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out;
  int mid, slow, size;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane, upto;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  size = nfast * sizeof(PACK_DATA);
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_plane;
    upto = slow * nmid * nfast;
    for (mid = 0; mid < nmid; mid++) {
      in = &(data[plane + mid * nstride_line]);
      out = &(buf[upto + mid * nfast]);
      memcpy(in, out, size);
    }
  }
}
static void unpack_3d_permute1_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane + mid]);
      end = begin + nfast * nstride_plane;
      for (in = begin; in < end; in += nstride_plane) *in = *(out++);
    }
  }
}
static void unpack_3d_permute1_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane + 2 * mid]);
      end = begin + nfast * nstride_plane;
      for (in = begin; in < end; in += nstride_plane) {
        *in = *(out++);
        *(in + 1) = *(out++);
      }
    }
  }
}
static void unpack_3d_permute1_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *instart, *begin, *end;
  int iqty, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, plane, nqty;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow * nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane + nqty * mid]);
      end = begin + nfast * nstride_plane;
      for (instart = begin; instart < end; instart += nstride_plane) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
      }
    }
  }
}
static void unpack_3d_permute2_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[slow + mid * nstride_plane]);
      end = begin + nfast * nstride_line;
      for (in = begin; in < end; in += nstride_line) *in = *(out++);
    }
  }
}
static void unpack_3d_permute2_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *begin, *end;
  int mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[2 * slow + mid * nstride_plane]);
      end = begin + nfast * nstride_line;
      for (in = begin; in < end; in += nstride_line) {
        *in = *(out++);
        *(in + 1) = *(out++);
      }
    }
  }
}
static void unpack_3d_permute2_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan)
{
  PACK_DATA *in, *out, *instart, *begin, *end;
  int iqty, mid, slow;
  int nfast, nmid, nslow, nstride_line, nstride_plane, nqty;
  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;
  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[nqty * slow + mid * nstride_plane]);
      end = begin + nfast * nstride_line;
      for (instart = begin; instart < end; instart += nstride_line) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
      }
    }
  }
}
#endif
