#include "error.h"
#include "input.h"
#include "universe.h"
using namespace LAMMPS_NS;
static std::string truncpath(const std::string &path) {
  std::size_t found = path.find("src/");
  if (found != std::string::npos)
    return path.substr(found);
  else
    return path;
}
Error::Error(LAMMPS *lmp)
    : Pointers(lmp), numwarn(0), maxwarn(100), allwarn(0) {}
void Error::universe_all(const std::string &file, int line,
                         const std::string &str) {
  MPI_Barrier(universe->uworld);
  std::string mesg = "ERROR: " + str;
  try {
    mesg += fmt::format(" ({}:{})\n", truncpath(file), line);
  } catch (fmt::format_error &) {
    ;
  }
  if (universe->me == 0) {
    if (universe->uscreen)
      fputs(mesg.c_str(), universe->uscreen);
    if (universe->ulogfile)
      fputs(mesg.c_str(), universe->ulogfile);
  }
  if (universe->nworlds > 1) {
    if (screen && screen != stdout)
      fclose(screen);
    if (logfile)
      fclose(logfile);
  }
  if (universe->ulogfile)
    fclose(universe->ulogfile);
  MPI_Finalize();
  exit(1);
}
void Error::universe_one(const std::string &file, int line,
                         const std::string &str) {
  std::string mesg = fmt::format("ERROR on proc {}: {} ({}:{})\n", universe->me,
                                 str, truncpath(file), line);
  if (universe->uscreen)
    fputs(mesg.c_str(), universe->uscreen);
  MPI_Abort(universe->uworld, 1);
  exit(1);
}
void Error::universe_warn(const std::string &file, int line,
                          const std::string &str) {
  ++numwarn;
  if ((maxwarn != 0) &&
      ((numwarn > maxwarn) || (allwarn > maxwarn) || (maxwarn < 0)))
    return;
  if (universe->uscreen)
    fmt::print(universe->uscreen, "WARNING on proc {}: {} ({}:{})\n",
               universe->me, str, truncpath(file), line);
}
void Error::all(const std::string &file, int line, const std::string &str) {
  MPI_Barrier(world);
  int me;
  std::string lastcmd = "(unknown)";
  MPI_Comm_rank(world, &me);
  if (me == 0) {
    std::string mesg = "ERROR: " + str;
    if (input && input->line)
      lastcmd = input->line;
    try {
      mesg += fmt::format(" ({}:{})\nLast command: {}\n", truncpath(file), line,
                          lastcmd);
    } catch (fmt::format_error &) {
      ;
    }
    utils::logmesg(lmp, mesg);
  }
  if (screen && screen != stdout)
    fclose(screen);
  if (logfile)
    fclose(logfile);
  if (universe->nworlds > 1)
    MPI_Abort(universe->uworld, 1);
  MPI_Finalize();
  exit(1);
}
void Error::one(const std::string &file, int line, const std::string &str) {
  int me;
  std::string lastcmd = "(unknown)";
  MPI_Comm_rank(world, &me);
  if (input && input->line)
    lastcmd = input->line;
  std::string mesg =
      fmt::format("ERROR on proc {}: {} ({}:{})\nLast command: {}\n", me, str,
                  truncpath(file), line, lastcmd);
  utils::logmesg(lmp, mesg);
  if (universe->nworlds > 1)
    if (universe->uscreen)
      fputs(mesg.c_str(), universe->uscreen);
  utils::flush_buffers(lmp);
  MPI_Abort(world, 1);
  exit(1);
}
void Error::_all(const std::string &file, int line, fmt::string_view format,
                 fmt::format_args args) {
  try {
    all(file, line, fmt::vformat(format, args));
  } catch (fmt::format_error &e) {
    all(file, line, e.what());
  }
  exit(1);
}
void Error::_one(const std::string &file, int line, fmt::string_view format,
                 fmt::format_args args) {
  try {
    one(file, line, fmt::vformat(format, args));
  } catch (fmt::format_error &e) {
    one(file, line, e.what());
  }
  exit(1);
}
void Error::warning(const std::string &file, int line, const std::string &str) {
  ++numwarn;
  if ((maxwarn != 0) &&
      ((numwarn > maxwarn) || (allwarn > maxwarn) || (maxwarn < 0)))
    return;
  std::string mesg =
      fmt::format("WARNING: {} ({}:{})\n", str, truncpath(file), line);
  if (screen)
    fputs(mesg.c_str(), screen);
  if (logfile)
    fputs(mesg.c_str(), logfile);
}
void Error::_warning(const std::string &file, int line, fmt::string_view format,
                     fmt::format_args args) {
  try {
    warning(file, line, fmt::vformat(format, args));
  } catch (fmt::format_error &e) {
    warning(file, line, e.what());
  }
}
void Error::message(const std::string &file, int line, const std::string &str) {
  std::string mesg = fmt::format("{} ({}:{})\n", str, truncpath(file), line);
  if (screen)
    fputs(mesg.c_str(), screen);
  if (logfile)
    fputs(mesg.c_str(), logfile);
}
void Error::_message(const std::string &file, int line, fmt::string_view format,
                     fmt::format_args args) {
  try {
    message(file, line, fmt::vformat(format, args));
  } catch (fmt::format_error &e) {
    message(file, line, e.what());
  }
}
void Error::done(int status) {
  MPI_Barrier(world);
  if (screen && screen != stdout)
    fclose(screen);
  if (logfile)
    fclose(logfile);
  MPI_Finalize();
  exit(status);
}
