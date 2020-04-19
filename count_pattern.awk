{
  count = 0
  $0 = tolower($0)
  pattern = tolower(pattern)
  while (length() > 0) {
    m = match($0, pattern)
    if (m == 0)
      break
    count++
    $0 = substr($0, m + 1)
  }
  print count
}