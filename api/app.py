import json
from http.server import BaseHTTPRequestHandler

def is_nucleotide(seq):
    nucleotides = set('ATGC')
    return all(nuc in nucleotides for nuc in seq)

def alignment_score(alignment1, alignment2, match_score, mismatch_score, gap_penalty):
    score = 0
    for char1, char2 in zip(alignment1, alignment2):
        if char1 == '-' or char2 == '-':
            score += gap_penalty
        elif char1 == char2:
            score += match_score
        else:
            score += mismatch_score
    return score

def needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty):

    matrix = [[0 for i in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]

    for i in range(1, len(seq1) + 1):
        matrix[i][0] = i * gap_penalty
    for j in range(1, len(seq2) + 1):
        matrix[0][j] = j * gap_penalty

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):

            a = matrix[i - 1][j] + gap_penalty
            b = matrix[i][j - 1] + gap_penalty
            c = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)

            matrix[i][j] = max(a, b, c)

    align1 = ''
    align2 = ''

    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:

        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1

        elif j > 0 and matrix[i][j] == matrix[i][j - 1] + gap_penalty:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

        else:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1

    return align1, align2, matrix


class handler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()

        response = {"message": "Needleman-Wunsch API running"}
        self.wfile.write(json.dumps(response).encode())
    
    def do_POST(self):

        content_length = int(self.headers['Content-Length'])
        body = self.rfile.read(content_length)

        data = json.loads(body)

        seq1 = data["seq1"].upper()
        seq2 = data["seq2"].upper()

        match_score = int(data["match"])
        mismatch_score = int(data["mismatch"])
        gap_penalty = int(data["gap"])

        alignment1, alignment2, matrix = needleman_wunsch(
            seq1, seq2, match_score, mismatch_score, gap_penalty
        )

        score = alignment_score(alignment1, alignment2, match_score, mismatch_score, gap_penalty)

        response = {
            "matrix": matrix,
            "alignment1": alignment1,
            "alignment2": alignment2,
            "score": score
        }

        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()

        self.wfile.write(json.dumps(response).encode())
