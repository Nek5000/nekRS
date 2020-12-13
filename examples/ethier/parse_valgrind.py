import re
def field_by_regex(regex,log_file_name, fieldnum = 0):
    with open(log_file_name) as log_file:
        content =log_file.readlines()
    content = [x.strip() for x in content]
    field=[]
    for line in content:
        m = re.search(regex,line)
        if m:
            field.append(int(m.groups()[fieldnum]))
    return field

def parse_valgrind():
  valgrind_regex = "ERROR SUMMARY: (.+?) errors from"
  logfile = "valgrind_output.log"
  errList = field_by_regex(valgrind_regex, logfile)
  assert len(errList) == 1
  numValgrindErrors = errList[0]
  if numValgrindErrors != 0:
    print(f"FAIL: found {numValgrindErrors} valgrind errors (supressed)")
    exit(2)
  else:
    print("PASS")
    exit(0)

if __name__ == "__main__":
  parse_valgrind()
